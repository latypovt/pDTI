#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Pig DTI BIDS pipeline (MRtrix + FSL bet4animal/dwifslpreproc)
#
# Usage examples:
#   pig_dti.sh --bids_dataset /path/to/bids --subject 9116
#   pig_dti.sh --bids_dataset /path/to/bids --subject 9116 --session Day14
#   pig_dti.sh --bids_dataset /path/to/bids --subject 9116 --deriv_name pigdti
#
# Assumes a BIDS-style layout like:
#   sub-9116/ses-Day14/anat/sub-9116_ses-Day14_run-01_T1w.nii.gz
#   sub-9116/ses-Day14/dwi/sub-9116_ses-Day14_dwi.nii.gz
#   sub-9116/ses-Day14/fmap/sub-9116_ses-Day14_epi.nii.gz
# ------------------------------------------------------------

show_help() {
  cat <<USAGE
Usage:
  $0 --bids_dataset PATH --subject ID [--session ID|all] [--deriv_name NAME]

Required:
  --bids_dataset PATH   Path to the root of the BIDS dataset
  --subject ID          Subject ID without 'sub-' prefix (e.g. 9116)

Optional:
  --session ID|all      Session label without 'ses-' prefix (e.g. Day14),
                        or 'all' to auto-detect all sessions with DWI.
                        Default: all
  --deriv_name NAME     Derivatives folder name (default: pigdti)

Expected per-session inputs:
  anat/
    sub-<ID>_ses-<SES>_run-01_T1w.nii[.gz]   (preferred)
    or sub-<ID>_ses-<SES>_T1w.nii[.gz]
    or any *T1w.nii[.gz] (as fallback)
  dwi/
    sub-<ID>_ses-<SES>_dwi.nii[.gz] + .bval/.bvec/.json
  fmap/ (optional, for reverse phase)
    sub-<ID>_ses-<SES>_epi.nii[.gz] + .json

Outputs:
  BIDS/derivatives/NAME/sub-<ID>/ses-<SES>/{anat,dwi,work}/...
  DTI metrics in T1 space:
    sub-<ID>_ses-<SES>_space-T1w_desc-FA/MD/AD/RD_dwi.nii.gz
USAGE
}

module load fsl
module load mrtrix3src
module load itksnap

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

# Prefer .nii, fall back to .nii.gz
choose_nii() {
  local base="$1"
  if [[ -f "${base}.nii" ]]; then
    echo "${base}.nii"
  elif [[ -f "${base}.nii.gz" ]]; then
    echo "${base}.nii.gz"
  else
    echo ""
  fi
}

# Choose a reasonable T1w from the anat folder
choose_t1w() {
  local anat_dir="$1"
  local sub_label="$2"   # e.g. sub-9116
  local ses_label="$3"   # e.g. ses-Day14

  local base_prefix="${anat_dir}/${sub_label}_${ses_label}"

  # 1) run-01_T1w
  local cand
  local nii

  cand="${base_prefix}_run-01_T1w"
  nii=$(choose_nii "${cand}")
  if [[ -n "${nii}" ]]; then
    echo "${nii}"
    return 0
  fi

  # 2) plain T1w (no run/acq)
  cand="${base_prefix}_T1w"
  nii=$(choose_nii "${cand}")
  if [[ -n "${nii}" ]]; then
    echo "${nii}"
    return 0
  fi

  # 3) specific mp2rage magnitude (if you prefer gw)
  cand="${base_prefix}_acq-mp2rageGw_T1w"
  nii=$(choose_nii "${cand}")
  if [[ -n "${nii}" ]]; then
    echo "${nii}"
    return 0
  fi

  # 4) fallback: first T1w in folder
  nii=$(ls "${anat_dir}"/${sub_label}_${ses_label}_*T1w.nii* 2>/dev/null | head -n1 || true)
  echo "${nii}"
}

# ------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------
BIDS_ROOT=""
SUB=""
SES="all"
DERIV_NAME="pigdti"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bids_dataset)
      BIDS_ROOT="$2"
      shift 2
      ;;
    --subject)
      SUB="$2"
      shift 2
      ;;
    --session)
      SES="$2"
      shift 2
      ;;
    --deriv_name)
      DERIV_NAME="$2"
      shift 2
      ;;
    -h|--help)
      show_help
      exit 0
      ;;
    *)
      echo "Unknown argument: $1"
      show_help
      exit 1
      ;;
  esac
done

if [[ -z "${BIDS_ROOT}" || -z "${SUB}" ]]; then
  echo "Error: --bids_dataset and --subject are required."
  show_help
  exit 1
fi

SUB_LABEL="sub-${SUB}"

if [[ ! -d "${BIDS_ROOT}/${SUB_LABEL}" ]]; then
  echo "ERROR: subject folder not found: ${BIDS_ROOT}/${SUB_LABEL}" >&2
  exit 1
fi

DERIV_ROOT="${BIDS_ROOT}/derivatives/${DERIV_NAME}"

# ------------------------------------------------------------
# Discover sessions (based on *_dwi.nii[.gz])
# ------------------------------------------------------------
SESSIONS=()

if [[ "${SES}" == "all" ]]; then
  echo "[INFO] Auto-detecting sessions for ${SUB_LABEL} with DWI ..."
  for ses_dir in "${BIDS_ROOT}/${SUB_LABEL}"/ses-*; do
    [[ -d "${ses_dir}" ]] || continue
    ses_base=$(basename "${ses_dir}")   # e.g. ses-Day14
    ses_id="${ses_base#ses-}"          # Day14

    dwi_dir="${ses_dir}/dwi"
    dwi_base="${dwi_dir}/${SUB_LABEL}_${ses_base}_dwi"
    dwi_nii=$(choose_nii "${dwi_base}")

    if [[ -n "${dwi_nii}" ]]; then
      echo "  + Found DWI for session ${ses_base}: $(basename "${dwi_nii}")"
      SESSIONS+=("${ses_id}")
    else
      echo "  - No DWI for session ${ses_base}, skipping"
    fi
  done
else
  SESSIONS+=("${SES}")
fi

if ((${#SESSIONS[@]} == 0)); then
  echo "No sessions with DWI found for ${SUB_LABEL}. Nothing to do."
  exit 0
fi

# ------------------------------------------------------------
# Session-level pipeline
# ------------------------------------------------------------
run_session() {
  local SES_LOCAL="$1"
  local SES_LABEL="ses-${SES_LOCAL}"

  echo "======================================================="
  echo "[INFO] Processing ${SUB_LABEL} ${SES_LABEL}"
  echo "======================================================="

  local SES_ROOT="${BIDS_ROOT}/${SUB_LABEL}/${SES_LABEL}"
  local ANAT_DIR="${SES_ROOT}/anat"
  local DWI_DIR="${SES_ROOT}/dwi"
  local FMAP_DIR="${SES_ROOT}/fmap"

  if [[ ! -d "${ANAT_DIR}" ]]; then
    echo "  [WARN] No anat folder for ${SUB_LABEL} ${SES_LABEL}, skipping."
    return
  fi
  if [[ ! -d "${DWI_DIR}" ]]; then
    echo "  [WARN] No dwi folder for ${SUB_LABEL} ${SES_LABEL}, skipping."
    return
  fi

  # ---- Choose T1w ----
  local T1W
  T1W=$(choose_t1w "${ANAT_DIR}" "${SUB_LABEL}" "${SES_LABEL}")

  if [[ -z "${T1W}" ]]; then
    echo "  [WARN] No T1w found in ${ANAT_DIR}, skipping session."
    return
  else
    echo "  [INFO] Using T1w image: $(basename "${T1W}")"
  fi

  # ---- DWI ----
  local DWI_BASE="${DWI_DIR}/${SUB_LABEL}_${SES_LABEL}_dwi"
  local DWI_NII
  DWI_NII=$(choose_nii "${DWI_BASE}")

  if [[ -z "${DWI_NII}" ]]; then
    echo "  [WARN] DWI not found for ${SUB_LABEL} ${SES_LABEL}, skipping."
    return
  fi

  local BVAL="${DWI_BASE}.bval"
  local BVEC="${DWI_BASE}.bvec"
  local DWI_JSON="${DWI_BASE}.json"

  if [[ ! -f "${BVAL}" || ! -f "${BVEC}" ]]; then
    echo "  [WARN] bval/bvec missing for DWI (${BVAL}, ${BVEC}), skipping."
    return
  fi

  # ---- fmap / epi (optional) ----
  local EPI_BASE="${FMAP_DIR}/${SUB_LABEL}_${SES_LABEL}_epi"
  local EPI_NII=""
  local EPI_JSON=""
  local HAVE_EPI="0"

  if [[ -d "${FMAP_DIR}" ]]; then
    EPI_NII=$(choose_nii "${EPI_BASE}")
    EPI_JSON="${EPI_BASE}.json"
    if [[ -n "${EPI_NII}" && -f "${EPI_JSON}" ]]; then
      HAVE_EPI="1"
      echo "  [INFO] Found SE-EPI fmap: $(basename "${EPI_NII}")"
    else
      echo "  [INFO] No usable EPI fmap; will run dwifslpreproc with -rpe_none."
    fi
  else
    echo "  [INFO] No fmap folder; will run dwifslpreproc with -rpe_none."
  fi

  local ANAT_DERIV="${DERIV_ROOT}/${SUB_LABEL}/${SES_LABEL}/anat"
  local WORK="${DERIV_ROOT}/${SUB_LABEL}/${SES_LABEL}/work"
  local DWI_DERIV="${DERIV_ROOT}/${SUB_LABEL}/${SES_LABEL}/dwi"

  mkdir -p "${ANAT_DERIV}" "${WORK}" "${DWI_DERIV}"

  # --------------------------------------------------------
  # 1) T1 brain mask with human-in-the-loop confirmation
  #    - Initial mask from bet4animal (-z 7 -R)
  #    - Launch ITK-SNAP for manual QC/editing
  #    - Apply morphological cleanup (fill holes, closing, smooth edges)
  #    - Only then create final skull-stripped T1 via fslmaths
  # --------------------------------------------------------
  local T1_BRAIN_NII="${ANAT_DERIV}/${SUB_LABEL}_${SES_LABEL}_T1w_brain.nii.gz"
  local T1_MASK_NII="${ANAT_DERIV}/${SUB_LABEL}_${SES_LABEL}_T1w_brain_mask.nii.gz"

  if [[ -f "${T1_BRAIN_NII}" && -f "${T1_MASK_NII}" ]]; then
    echo "  [INFO] T1w brain and mask already exist; assuming previously confirmed. Skipping ITK-SNAP + morphology."
  else
    # 1a) Generate initial mask if missing
    if [[ ! -f "${T1_MASK_NII}" ]]; then
      echo "  [INFO] Running bet4animal for initial brain mask (-z 7 -R) ..."
      local T1_BRAIN_INIT="${ANAT_DERIV}/${SUB_LABEL}_${SES_LABEL}_T1w_brain_init.nii.gz"
      bet4animal "${T1W}" "${T1_BRAIN_INIT}" -m -z 7 -R

      local T1_MASK_INIT="${ANAT_DERIV}/${SUB_LABEL}_${SES_LABEL}_T1w_brain_init_mask.nii.gz"
      if [[ ! -f "${T1_MASK_INIT}" ]]; then
        echo "  [ERROR] bet4animal did not produce expected mask: ${T1_MASK_INIT}" >&2
        exit 1
      fi

      # Use this as the editable mask
      cp "${T1_MASK_INIT}" "${T1_MASK_NII}"
    else
      echo "  [INFO] Found existing T1 mask: $(basename "${T1_MASK_NII}")"
    fi

    # 1b) Check ITK-SNAP availability
    if ! command -v itksnap >/dev/null 2>&1; then
      echo "  [ERROR] itksnap not found in PATH; cannot perform interactive mask editing." >&2
      echo "         Either install ITK-SNAP or precompute ${T1_MASK_NII} and rerun." >&2
      exit 1
    fi

    # 1c) Launch ITK-SNAP for interactive editing
    echo "  [INFO] Launching ITK-SNAP for brain mask QC/editing."
    echo "         Grey image     : ${T1W}"
    echo "         Segmentation   : ${T1_MASK_NII}"
    echo "         Please:"
    echo "           - Inspect and edit the mask as needed."
    echo "           - SAVE the segmentation (overwrite ${T1_MASK_NII})."
    echo "           - Close ITK-SNAP to continue the pipeline."
    itksnap -g "${T1W}" -s "${T1_MASK_NII}"

    # 1d) Morphological cleanup: fill holes, closing, smooth edges
    echo "  [INFO] Post-edit mask cleanup: fill holes, closing, edge smoothing ..."
    local T1_MASK_TMP="${ANAT_DERIV}/${SUB_LABEL}_${SES_LABEL}_T1w_brain_mask_tmp.nii.gz"

    # Fill internal holes (3D)
    fslmaths "${T1_MASK_NII}" -fillh "${T1_MASK_TMP}"

    # Binary closing: dilate then erode to close small gaps
    fslmaths "${T1_MASK_TMP}" -dilM -ero "${T1_MASK_TMP}"

    # Smooth edges slightly, then re-binarize
    # (sigma=0.5 is gentle; tweak if needed)
    fslmaths "${T1_MASK_TMP}" -s 0.5 -thr 0.5 -bin "${T1_MASK_NII}"

    # 1e) Create final skull-stripped T1 from cleaned mask
    echo "  [INFO] Creating final skull-stripped T1 from cleaned mask ..."
    fslmaths "${T1W}" -mas "${T1_MASK_NII}" "${T1_BRAIN_NII}"
  fi

  # --------------------------------------------------------
  # 2) Convert DWI → .mif (skip if exists)
  # --------------------------------------------------------
  if [[ ! -f "${WORK}/dwi.mif" ]]; then
    echo "  [INFO] Converting DWI to MRtrix format ..."
    if [[ -f "${DWI_JSON}" ]]; then
      mrconvert "${DWI_NII}" "${WORK}/dwi.mif" \
        -fslgrad "${BVEC}" "${BVAL}" \
        -json_import "${DWI_JSON}" \
        -datatype float32 -stride 0,0,0,1
    else
      echo "  [WARN] DWI JSON not found (${DWI_JSON}), importing without JSON."
      mrconvert "${DWI_NII}" "${WORK}/dwi.mif" \
        -fslgrad "${BVEC}" "${BVAL}" \
        -datatype float32 -stride 0,0,0,1
    fi
  else
    echo "  [INFO] dwi.mif already exists, skipping mrconvert."
  fi

  # --------------------------------------------------------
  # 3) Denoising + Gibbs (skip if outputs exist)
  # --------------------------------------------------------
  if [[ ! -f "${WORK}/dwi_den_degibbs.mif" ]]; then
    echo "  [INFO] Denoising (dwidenoise) ..."
    dwidenoise "${WORK}/dwi.mif" "${WORK}/dwi_den.mif" -noise "${WORK}/noise.mif"

    echo "  [INFO] Gibbs ringing removal (mrdegibbs) ..."
    mrdegibbs "${WORK}/dwi_den.mif" "${WORK}/dwi_den_degibbs.mif"
  else
    echo "  [INFO] Denoising + Gibbs already done, skipping."
  fi

  # --------------------------------------------------------
  # 4) TOPUP / Eddy (dwifslpreproc)
  #    If SE-EPI fmap exists -> use -rpe_pair; else -> -rpe_none.
  # --------------------------------------------------------
  if [[ ! -f "${WORK}/dwi_preproc.mif" ]]; then
    if [[ "${HAVE_EPI}" == "1" ]]; then
      echo "  [INFO] Converting SE-EPI fmap to .mif ..."
      if [[ ! -f "${WORK}/se_epi.mif" ]]; then
        mrconvert "${EPI_NII}" "${WORK}/se_epi.mif" \
          -json_import "${EPI_JSON}" \
          -datatype float32
      fi

      echo "  [INFO] Running dwifslpreproc with -rpe_pair and SE-EPI ..."
      dwifslpreproc "${WORK}/dwi_den_degibbs.mif" "${WORK}/dwi_preproc.mif" \
        -pe_dir ap \
        -rpe_pair \
        -se_epi "${WORK}/se_epi.mif" \
        -eddy_options " --slm=linear " \
        -nocleanup
    else
      echo "  [INFO] Running dwifslpreproc with -rpe_none ..."
      dwifslpreproc "${WORK}/dwi_den_degibbs.mif" "${WORK}/dwi_preproc.mif" \
        -pe_dir ap \
        -rpe_none \
        -eddy_options " --slm=linear " \
        -nocleanup
    fi
  else
    echo "  [INFO] dwi_preproc.mif exists, skipping dwifslpreproc."
  fi

  # --------------------------------------------------------
  # 5) Bias field correction (skip if biascorr exists)
  # --------------------------------------------------------
  if [[ ! -f "${WORK}/dwi_preproc_biascorr.mif" ]]; then
    echo "  [INFO] Bias field correction (dwibiascorrect ants) ..."
    dwibiascorrect ants \
      "${WORK}/dwi_preproc.mif" \
      "${WORK}/dwi_preproc_biascorr.mif" \
      -bias "${WORK}/bias.mif"
  else
    echo "  [INFO] dwi_preproc_biascorr.mif exists, skipping dwibiascorrect."
  fi

  # --------------------------------------------------------
  # 6) Convert T1 brain + mask to .mif (skip if exists)
  # --------------------------------------------------------
  if [[ ! -f "${WORK}/T1_brain.mif" ]]; then
    mrconvert "${T1_BRAIN_NII}" "${WORK}/T1_brain.mif"
  fi
  if [[ ! -f "${WORK}/T1_mask.mif" ]]; then
    mrconvert "${T1_MASK_NII}"  "${WORK}/T1_mask.mif"
  fi

  # --------------------------------------------------------
  # 7) Resample corrected DWI to T1 grid (NO registration)
  # --------------------------------------------------------
  if [[ ! -f "${WORK}/dwi_inT1.mif" ]]; then
    echo "  [INFO] Resampling corrected DWI to T1w resolution/grid ..."
    mrgrid "${WORK}/dwi_preproc_biascorr.mif" regrid \
      "${WORK}/dwi_inT1.mif" \
      -template "${WORK}/T1_brain.mif"
  else
    echo "  [INFO] dwi_inT1.mif exists, skipping mrgrid."
  fi

  # --------------------------------------------------------
  # 8) Tensor fitting + metrics in T1 space using T1 mask
  # --------------------------------------------------------
  if [[ ! -f "${WORK}/dt_T1.mif" ]]; then
    echo "  [INFO] Tensor fitting (dwi2tensor) in T1 space using T1 mask ..."
    dwi2tensor "${WORK}/dwi_inT1.mif" "${WORK}/dt_T1.mif" \
      -mask "${WORK}/T1_mask.mif"
  else
    echo "  [INFO] dt_T1.mif exists, skipping dwi2tensor."
  fi

  if [[ ! -f "${WORK}/fa_T1.mif" ]]; then
    echo "  [INFO] Computing tensor metrics (FA/MD/AD/RD) in T1 space ..."
    tensor2metric "${WORK}/dt_T1.mif" \
      -fa  "${WORK}/fa_T1.mif" \
      -adc "${WORK}/md_T1.mif" \
      -ad  "${WORK}/ad_T1.mif" \
      -rd  "${WORK}/rd_T1.mif"
  else
    echo "  [INFO] T1-space tensor metrics already exist, skipping tensor2metric."
  fi

  # --------------------------------------------------------
  # 9) Export NIfTI to derivatives (T1-space metrics)
  # --------------------------------------------------------
  if [[ ! -f "${DWI_DERIV}/${SUB_LABEL}_${SES_LABEL}_space-T1w_desc-FA_dwi.nii.gz" ]]; then
    echo "  [INFO] Exporting T1w-space DTI scalars ..."
    mrconvert "${WORK}/fa_T1.mif" \
      "${DWI_DERIV}/${SUB_LABEL}_${SES_LABEL}_space-T1w_desc-FA_dwi.nii.gz"
    mrconvert "${WORK}/md_T1.mif" \
      "${DWI_DERIV}/${SUB_LABEL}_${SES_LABEL}_space-T1w_desc-MD_dwi.nii.gz"
    mrconvert "${WORK}/ad_T1.mif" \
      "${DWI_DERIV}/${SUB_LABEL}_${SES_LABEL}_space-T1w_desc-AD_dwi.nii.gz"
    mrconvert "${WORK}/rd_T1.mif" \
      "${DWI_DERIV}/${SUB_LABEL}_${SES_LABEL}_space-T1w_desc-RD_dwi.nii.gz"

    mrconvert "${WORK}/T1_mask.mif" \
      "${DWI_DERIV}/${SUB_LABEL}_${SES_LABEL}_space-T1w_desc-brain_mask.nii.gz"
  else
    echo "  [INFO] T1-space DTI NIfTIs already exist, skipping export."
  fi

  echo "  [INFO] Finished ${SUB_LABEL} ${SES_LABEL}"
}

# ------------------------------------------------------------
# Loop over sessions
# ------------------------------------------------------------
for ses_id in "${SESSIONS[@]}"; do
  run_session "${ses_id}"
done

echo "[INFO] All done for ${SUB_LABEL}. Derivatives in:"
echo "  ${DERIV_ROOT}/${SUB_LABEL}/ses-*/anat"
echo "  ${DERIV_ROOT}/${SUB_LABEL}/ses-*/dwi"
