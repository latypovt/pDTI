#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Pig DTI BIDS pipeline: Updated with JSON Readout & Mask Polishing
# ------------------------------------------------------------

# --- MODULES ---
module load fsl 2>/dev/null || true
module load mrtrix3src 2>/dev/null || true
module load itksnap 2>/dev/null || true

# --- PIG REORIENTATION ---
PERMUTE_AXES="0,2,1"
FLIP_AXES="2"

# --- ARG PARSING ---
BIDS_ROOT=""; SUB=""; SES="all"; DERIV_NAME="pdti2"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bids_dataset) BIDS_ROOT="$2"; shift 2 ;;
    --subject) SUB="$2"; shift 2 ;;
    --session) SES="$2"; shift 2 ;;
    --deriv_name) DERIV_NAME="$2"; shift 2 ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

SUB_LABEL="sub-${SUB}"
DERIV_ROOT="${BIDS_ROOT}/derivatives/${DERIV_NAME}"
SESSIONS=($(ls -d ${BIDS_ROOT}/${SUB_LABEL}/ses-* 2>/dev/null | xargs -n1 basename | sed 's/ses-//'))
[[ "${SES}" != "all" ]] && SESSIONS=("${SES}")

run_session() {
  local SES_LABEL="ses-$1"
  local WORK="${DERIV_ROOT}/${SUB_LABEL}/${SES_LABEL}/work"
  local DWI_DERIV="${DERIV_ROOT}/${SUB_LABEL}/${SES_LABEL}/dwi"
  local ANAT_DERIV="${DERIV_ROOT}/${SUB_LABEL}/${SES_LABEL}/anat"
  mkdir -p "${WORK}" "${DWI_DERIV}" "${ANAT_DERIV}"

  local T1W="${BIDS_ROOT}/${SUB_LABEL}/${SES_LABEL}/anat/${SUB_LABEL}_${SES_LABEL}_run-01_T1w.nii.gz"
  local DWI_BASE="${BIDS_ROOT}/${SUB_LABEL}/${SES_LABEL}/dwi/${SUB_LABEL}_${SES_LABEL}_dwi"
  local DWI_JSON="${BIDS_ROOT}/${SUB_LABEL}/${SES_LABEL}/dwi/${SUB_LABEL}_${SES_LABEL}_dwi.json"

  # 1) PREPROCESSING DWI (Raw Space)
  if [[ ! -f "${WORK}/1_dwi_preproc.mif" ]]; then
    echo "  [CHECKPOINT 1] Preprocessing in scanner space..."
    READOUT=$(jq -r '.EstimatedTotalReadoutTime' "${DWI_JSON}")
    mrconvert "${DWI_BASE}.nii.gz" -fslgrad "${DWI_BASE}.bvec" "${DWI_BASE}.bval" "${WORK}/raw.mif" -force
    dwidenoise "${WORK}/raw.mif" "${WORK}/den.mif" -force
    mrdegibbs "${WORK}/den.mif" "${WORK}/degibbs.mif" -force
    dwifslpreproc "${WORK}/degibbs.mif" "${WORK}/eddy.mif" -rpe_none -pe_dir ap \
      -readout_time "${READOUT}" -eddy_options " --slm=linear " -force
    dwibiascorrect ants "${WORK}/eddy.mif" "${WORK}/1_dwi_preproc.mif" -force
  fi

  # 2) RESAMPLING DWI TO T1 RESOLUTION
  if [[ ! -f "${WORK}/2_dwi_resampled.mif" ]]; then
    echo "  [CHECKPOINT 2] Resampling DWI to T1 grid..."
    mrgrid "${WORK}/1_dwi_preproc.mif" regrid -template "${T1W}" -interp sinc "${WORK}/2_dwi_resampled.mif" -force
  fi

  # 3) FLIPPING AND PERMUTING AXES
  if [[ ! -f "${WORK}/3_T1_reo.nii.gz" || ! -f "${WORK}/3_dwi_reo.mif" ]]; then
    echo "  [CHECKPOINT 3] Reorienting to Pig-on-Belly space..."
    T1_DIMS=$(mrinfo "${T1W}" -ndim)
    if [ "$T1_DIMS" -eq 4 ]; then
      mrconvert "${T1W}" -coord 3 0 -axes "${PERMUTE_AXES}" - | mrtransform - "${WORK}/3_T1_reo.nii.gz" -flip "${FLIP_AXES}" -force
    else
      mrconvert "${T1W}" -axes "${PERMUTE_AXES}" - | mrtransform - "${WORK}/3_T1_reo.nii.gz" -flip "${FLIP_AXES}" -force
    fi
    mrconvert "${WORK}/2_dwi_resampled.mif" -axes "${PERMUTE_AXES},3" - | mrtransform - "${WORK}/3_dwi_reo.mif" -flip "${FLIP_AXES}" -force
  fi

  # 4) DUAL BRAIN MASKING & POLISHING
  local T1_MASK="${WORK}/4_T1_mask.nii.gz"
  local DWI_MASK="${WORK}/4_DWI_mask.nii.gz"
  
  if [[ ! -f "${T1_MASK}" || ! -f "${DWI_MASK}" ]]; then
    echo "  [CHECKPOINT 4] Generating and Polishing Masks..."
    
    # Initial Mask Generation
    bet4animal "${WORK}/3_T1_reo.nii.gz" "${WORK}/T1_brain_init.nii.gz" -m -z 7 -R
    mv "${WORK}/T1_brain_init_mask.nii.gz" "${T1_MASK}"
    
    dwiextract "${WORK}/3_dwi_reo.mif" - -bzero | mrmath - mean "${WORK}/3_meanb0_reo.nii.gz" -axis 3 -force
    bet4animal "${WORK}/3_meanb0_reo.nii.gz" "${WORK}/DWI_brain_init.nii.gz" -m -z 7 -R
    mv "${WORK}/DWI_brain_init_mask.nii.gz" "${DWI_MASK}"

    echo "  [USER ACTION] Edit masks in ITK-SNAP..."
    itksnap -g "${WORK}/3_T1_reo.nii.gz" -s "${T1_MASK}"
    itksnap -g "${WORK}/3_meanb0_reo.nii.gz" -s "${DWI_MASK}"
    
    # Polishing Step: Largest Connected Component + Edge Smoothing
    for MASK in "${T1_MASK}" "${DWI_MASK}"; do
      echo "    Polishing ${MASK}..."
      # 1. Keep only the largest connected component (removes isolated <10 voxel blobs)
      maskfilter "${MASK}" connect -largest "${WORK}/tmp_connected.nii.gz" -force
      # 2. Smooth edges (Median filter to "polish" borders without losing volume)
      maskfilter "${WORK}/tmp_connected.nii.gz" median "${MASK}" -force
      # 3. Final binary check and hole filling
      fslmaths "${MASK}" -fillh -bin "${MASK}"
    done

    # Export to derivatives
    cp "${WORK}/3_T1_reo.nii.gz" "${ANAT_DERIV}/${SUB_LABEL}_${SES_LABEL}_desc-preproc_T1w.nii.gz"
    fslmaths "${ANAT_DERIV}/${SUB_LABEL}_${SES_LABEL}_desc-preproc_T1w.nii.gz" -mas "${T1_MASK}" "${ANAT_DERIV}/${SUB_LABEL}_${SES_LABEL}_desc-brain_T1w.nii.gz"
  fi

  # 5) DTI FITTING
  if [[ ! -f "${WORK}/5_dt.mif" ]]; then
    echo "  [CHECKPOINT 5] Fitting tensors using manual polished DWI mask..."
    dwi2tensor "${WORK}/3_dwi_reo.mif" "${WORK}/5_dt.mif" -mask "${DWI_MASK}" -force
  fi

  # 6) SCALARS EXPORT
  echo "  [CHECKPOINT 6] Exporting scalars..."
  tensor2metric "${WORK}/5_dt.mif" \
    -fa "${DWI_DERIV}/${SUB_LABEL}_${SES_LABEL}_desc-FA_dwi.nii.gz" \
    -adc "${DWI_DERIV}/${SUB_LABEL}_${SES_LABEL}_desc-MD_dwi.nii.gz" \
    -ad "${DWI_DERIV}/${SUB_LABEL}_${SES_LABEL}_desc-AD_dwi.nii.gz" \
    -rd "${DWI_DERIV}/${SUB_LABEL}_${SES_LABEL}_desc-RD_dwi.nii.gz" -force

  echo "  [SUCCESS] Subject ${SUB_LABEL} sequence complete."
}

for s in "${SESSIONS[@]}"; do run_session "$s"; done