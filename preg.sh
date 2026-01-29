#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Pig Atlas Registration Pipeline (Non-linear)
# 1. T1_brain (subject) -> T1_template (atlas)
# 2. FA (subject) -> FA_template (atlas) + apply to MD/AD/RD
# ------------------------------------------------------------

log(){ echo "[$(date +'%F %T')] $*"; }
die(){ echo "ERROR: $*" >&2; exit 1; }
need(){ command -v "$1" >/dev/null 2>&1 || die "Missing command: $1"; }

usage() {
  cat <<EOF
Usage:
  $(basename "$0") --bids_dataset PATH --subject ID [options]

Required:
  --bids_dataset PATH    Path to BIDS root
  --subject ID           Subject ID (without sub-)

Options:
  --session SES          Process specific session (default: all)
  --deriv_name NAME      Input DTI derivative name (default: pdti2)
  --out_deriv NAME       Output derivative name (default: atlasreg)
  --atlas_dir PATH       Path to atlas folder (default: atlas/4wkAtlas)
  --transform_type TYPE  ANTs registration type (default: SyN)
  --python_venv PATH     Path to python executable with ANTs installed
EOF
}

# Default Paths
DERIV_NAME="pdti2"
OUT_DERIV="atlasreg"
TRANSFORM_TYPE="SyN"
ATLAS_DIR="atlas/4wkAtlas"
PYTHON_VENV="/home/timurlatypov/.virtualenvs/tim_apple/bin/python" # User should override if venv is needed

BIDS=""; SUB=""; SES=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bids_dataset) BIDS="$2"; shift 2 ;;
    --subject) SUB="$2"; shift 2 ;;
    --session) SES="$2"; shift 2 ;;
    --deriv_name) DERIV_NAME="$2"; shift 2 ;;
    --out_deriv) OUT_DERIV="$2"; shift 2 ;;
    --atlas_dir) ATLAS_DIR="$2"; shift 2 ;;
    --transform_type) TRANSFORM_TYPE="$2"; shift 2 ;;
    --python_venv) PYTHON_VENV="$2"; shift 2 ;;
    *) die "Unknown arg: $1" ;;
  esac
done

[[ -n "$BIDS" && -n "$SUB" ]] || { usage; die "Missing required args"; }

# Define Templates
T1_TEMPLATE="${ATLAS_DIR}/Pig_Brain_Atlas_4wk.nii"
FA_TEMPLATE="${ATLAS_DIR}/Pig_FA_Standard_4wk.nii"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANTSREG_PY="$SCRIPT_DIR/ANTsReg.py"

# Session discovery
SUB_DIR="$BIDS/sub-${SUB}"
if [[ -n "$SES" ]]; then
  SESSIONS=("ses-${SES}")
else
  mapfile -t SESSIONS < <(find "$SUB_DIR" -maxdepth 1 -type d -name "ses-*" -printf "%f\n" | sort)
fi

for sesdir in "${SESSIONS[@]}"; do
  log "---- Processing sub-${SUB}/${sesdir} ----"
  
  # Input directories from pdti2
  IN_ANAT_DIR="$BIDS/derivatives/${DERIV_NAME}/sub-${SUB}/${sesdir}/anat"
  IN_DWI_DIR="$BIDS/derivatives/${DERIV_NAME}/sub-${SUB}/${sesdir}/dwi"
  
  # Output directories
  WORK_DIR="$BIDS/derivatives/${OUT_DERIV}/sub-${SUB}/${sesdir}/work"
  OUT_ANAT_DIR="$BIDS/derivatives/${OUT_DERIV}/sub-${SUB}/${sesdir}/anat"
  OUT_DWI_DIR="$BIDS/derivatives/${OUT_DERIV}/sub-${SUB}/${sesdir}/dwi"
  mkdir -p "$WORK_DIR" "$OUT_ANAT_DIR" "$OUT_DWI_DIR"

  # Locate subject images
  shopt -s nullglob
  T1_BRAIN=("$IN_ANAT_DIR"/*_desc-brain_T1w.nii*)
  FA=("$IN_DWI_DIR"/*_desc-FA_dwi.nii*)
  MD=("$IN_DWI_DIR"/*_desc-MD_dwi.nii*)
  AD=("$IN_DWI_DIR"/*_desc-AD_dwi.nii*)
  RD=("$IN_DWI_DIR"/*_desc-RD_dwi.nii*)
  shopt -u nullglob

  # --- 1. T1 Registration (Subject T1_brain -> Atlas T1) ---
  if [[ -f "${T1_BRAIN[0]:-}" ]]; then
    log "Registering T1_brain to Atlas..."
    $PYTHON_VENV "$ANTSREG_PY" \
      --fixed_image "$T1_TEMPLATE" \
      --moving_image "${T1_BRAIN[0]}" \
      --output_prefix "${WORK_DIR}/t1_to_atlas" \
      --transform_type "$TRANSFORM_TYPE" \
      --apply_out_dir "$OUT_ANAT_DIR" \
      --apply_suffix "_space-Atlas" \
      --apply_images "${T1_BRAIN[0]}"
    
    # BIDS Rename
    mv "${OUT_ANAT_DIR}/$(basename "${T1_BRAIN[0]}" .nii.gz)_space-Atlas.nii.gz" \
       "${OUT_ANAT_DIR}/sub-${SUB}_${sesdir}_space-Atlas_desc-brain_T1w.nii.gz"
  fi

  # --- 2. FA Registration (Subject FA -> Atlas FA) ---
  if [[ -f "${FA[0]:-}" ]]; then
    log "Registering FA to Atlas + applying to metrics..."
    
    # Prepare list of scalar images to transform
    APPLY_LIST=()
    [[ -f "${MD[0]:-}" ]] && APPLY_LIST+=("${MD[0]}")
    [[ -f "${AD[0]:-}" ]] && APPLY_LIST+=("${AD[0]}")
    [[ -f "${RD[0]:-}" ]] && APPLY_LIST+=("${RD[0]}")

    $PYTHON_VENV "$ANTSREG_PY" \
      --fixed_image "$FA_TEMPLATE" \
      --moving_image "${FA[0]}" \
      --output_prefix "${WORK_DIR}/fa_to_atlas" \
      --transform_type "$TRANSFORM_TYPE" \
      --apply_out_dir "$OUT_DWI_DIR" \
      --apply_suffix "_space-Atlas" \
      --apply_images "${FA[0]}" "${APPLY_LIST[@]}"

    # BIDS Rename FA
    mv "${OUT_DWI_DIR}/$(basename "${FA[0]}" .nii.gz)_space-Atlas.nii.gz" \
       "${OUT_DWI_DIR}/sub-${SUB}_${sesdir}_space-Atlas_desc-FA_dwi.nii.gz"

    # BIDS Rename Scalars
    for metric_path in "${APPLY_LIST[@]}"; do
      m_name=$(basename "$metric_path" .nii.gz)
      metric_tag=$(echo "$m_name" | grep -oP "desc-\K[A-Z]{2}")
      mv "${OUT_DWI_DIR}/${m_name}_space-Atlas.nii.gz" \
         "${OUT_DWI_DIR}/sub-${SUB}_${sesdir}_space-Atlas_desc-${metric_tag}_dwi.nii.gz"
    done
  fi

  log "Done sub-${SUB}/${sesdir}"
done

log "Registration Pipeline Complete."