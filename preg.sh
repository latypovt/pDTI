#!/usr/bin/env bash
set -euo pipefail

log(){ echo "[$(date +'%F %T')] $*"; }
die(){ echo "ERROR: $*" >&2; exit 1; }
need(){ command -v "$1" >/dev/null 2>&1 || die "Missing command: $1"; }

usage() {
  cat <<EOF
Usage:
  $(basename "$0") --bids_dataset PATH --subject ID --atlas_fa PATH [--session SES] [options]

Required:
  --bids_dataset PATH
  --subject ID
  --atlas_fa PATH

Optional:
  --session SES                    If omitted, processes all ses-* for subject
  --deriv_name NAME                Input derivative (default: pigdti)
  --out_deriv NAME                 Output derivative (default: atlas_space)
  --permute_axes "0,2,1"           Default: 0,2,1
  --flip_axes "2"                  Default: 2
  --transform_type TYPE            Default: QuickRigid
  --dry_run
EOF
}

DERIV_NAME="pigdti"
OUT_DERIV="atlas_space"
PERMUTE_AXES="0,2,1"
FLIP_AXES="2"
TRANSFORM_TYPE="QuickRigid"
DRY_RUN=0

BIDS=""; SUB=""; SES=""; ATLAS_FA=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bids_dataset) BIDS="$2"; shift 2 ;;
    --subject) SUB="$2"; shift 2 ;;
    --session) SES="$2"; shift 2 ;;
    --atlas_fa) ATLAS_FA="$2"; shift 2 ;;
    --deriv_name) DERIV_NAME="$2"; shift 2 ;;
    --out_deriv) OUT_DERIV="$2"; shift 2 ;;
    --permute_axes) PERMUTE_AXES="$2"; shift 2 ;;
    --flip_axes) FLIP_AXES="$2"; shift 2 ;;
    --transform_type) TRANSFORM_TYPE="$2"; shift 2 ;;
    --dry_run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown arg: $1" ;;
  esac
done

[[ -n "$BIDS" && -n "$SUB" && -n "$ATLAS_FA" ]] || { usage; die "Missing required args"; }
[[ -d "$BIDS" ]] || die "No BIDS dir: $BIDS"
[[ -f "$ATLAS_FA" ]] || die "No atlas FA: $ATLAS_FA"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REORIENT_PY="$SCRIPT_DIR/reorient_image.py"
ANTSREG_PY="$SCRIPT_DIR/ANTsReg.py"
[[ -f "$REORIENT_PY" ]] || die "Missing $REORIENT_PY"
[[ -f "$ANTSREG_PY" ]] || die "Missing $ANTSREG_PY"

need python
need antsApplyTransforms

run() {
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY_RUN] $*"
  else
    eval "$@"
  fi
}

reoriented_path() {
  local in="$1"
  if [[ "$in" == *.nii.gz ]]; then
    echo "${in%.nii.gz}_reoriented.nii.gz"
  elif [[ "$in" == *.nii ]]; then
    echo "${in%.nii}_reoriented.nii"
  else
    die "Unsupported NIfTI extension: $in"
  fi
}

# Session discovery
SUB_DIR="$BIDS/sub-${SUB}"
[[ -d "$SUB_DIR" ]] || die "No subject dir: $SUB_DIR"

if [[ -n "$SES" ]]; then
  SESSIONS=("ses-${SES}")
else
  mapfile -t SESSIONS < <(find "$SUB_DIR" -maxdepth 1 -type d -name "ses-*" -printf "%f\n" | sort)
  [[ "${#SESSIONS[@]}" -gt 0 ]] || die "No ses-* found for sub-${SUB}"
fi

log "sub-${SUB} sessions=${SESSIONS[*]}"
log "permute_axes=$PERMUTE_AXES flip_axes=$FLIP_AXES"

for sesdir in "${SESSIONS[@]}"; do
  log "---- Processing sub-${SUB}/${sesdir} ----"

  IN_DWI_DIR="$BIDS/derivatives/${DERIV_NAME}/sub-${SUB}/${sesdir}/dwi"
  [[ -d "$IN_DWI_DIR" ]] || { log "Skip: missing $IN_DWI_DIR"; continue; }

  shopt -s nullglob
  FA_CAND=("$IN_DWI_DIR"/sub-"$SUB"_"$sesdir"_*desc-FA*_dwi.nii* "$IN_DWI_DIR"/sub-"$SUB"_"$sesdir"_*desc-FA*.nii*)
  MD_CAND=("$IN_DWI_DIR"/sub-"$SUB"_"$sesdir"_*desc-MD*_dwi.nii* "$IN_DWI_DIR"/sub-"$SUB"_"$sesdir"_*desc-MD*.nii*)
  AD_CAND=("$IN_DWI_DIR"/sub-"$SUB"_"$sesdir"_*desc-AD*_dwi.nii* "$IN_DWI_DIR"/sub-"$SUB"_"$sesdir"_*desc-AD*.nii*)
  RD_CAND=("$IN_DWI_DIR"/sub-"$SUB"_"$sesdir"_*desc-RD*_dwi.nii* "$IN_DWI_DIR"/sub-"$SUB"_"$sesdir"_*desc-RD*.nii*)
  shopt -u nullglob

  FA="${FA_CAND[0]:-}"
  MD="${MD_CAND[0]:-}"
  AD="${AD_CAND[0]:-}"
  RD="${RD_CAND[0]:-}"

  [[ -f "$FA" ]] || { log "Skip: FA not found"; continue; }

  # Reorient once (save into pigdti derivative with _reoriented suffix)
  FA_REO="$(reoriented_path "$FA")"
  MD_REO="$( [[ -f "${MD:-}" ]] && reoriented_path "$MD" || echo "" )"
  AD_REO="$( [[ -f "${AD:-}" ]] && reoriented_path "$AD" || echo "" )"
  RD_REO="$( [[ -f "${RD:-}" ]] && reoriented_path "$RD" || echo "" )"

  if [[ ! -f "$FA_REO" ]]; then
    log "Reorient FA -> $FA_REO"
    run "python \"$REORIENT_PY\" --in \"$FA\" --out \"$FA_REO\" --permute_axes \"$PERMUTE_AXES\" --flip_axes \"$FLIP_AXES\""
  else
    log "FA already reoriented: $FA_REO"
  fi

  for metric in MD AD RD; do
    src_var="${metric}"
    reo_var="${metric}_REO"
    src="${!src_var:-}"
    reo="${!reo_var:-}"
    [[ -n "$src" && -n "$reo" && -f "$src" ]] || { log "Skip $metric: missing"; continue; }
    if [[ ! -f "$reo" ]]; then
      log "Reorient $metric -> $reo"
      run "python \"$REORIENT_PY\" --in \"$src\" --out \"$reo\" --permute_axes \"$PERMUTE_AXES\" --flip_axes \"$FLIP_AXES\""
    else
      log "$metric already reoriented: $reo"
    fi
  done

  # Output derivatives/atlas_space
  OUT_DWI_DIR="$BIDS/derivatives/${OUT_DERIV}/sub-${SUB}/${sesdir}/dwi"
  run "mkdir -p \"$OUT_DWI_DIR\""

  OUT_PREFIX="$OUT_DWI_DIR/sub-${SUB}_${sesdir}_space-Atlas_desc-FA"

  # Register on reoriented FA and apply transforms immediately to MD/AD/RD (also reoriented)
  APPLY_LIST=()
  [[ -f "${MD_REO:-}" ]] && APPLY_LIST+=("\"$MD_REO\"")
  [[ -f "${AD_REO:-}" ]] && APPLY_LIST+=("\"$AD_REO\"")
  [[ -f "${RD_REO:-}" ]] && APPLY_LIST+=("\"$RD_REO\"")

  log "Register FA (reoriented) + apply to MD/AD/RD"
  if [[ "${#APPLY_LIST[@]}" -gt 0 ]]; then
    run "python \"$ANTSREG_PY\" \
      --fixed_image \"$ATLAS_FA\" \
      --moving_image \"$FA_REO\" \
      --output_prefix \"$OUT_PREFIX\" \
      --transform_type \"$TRANSFORM_TYPE\" \
      --apply_out_dir \"$OUT_DWI_DIR\" \
      --apply_suffix \"\" \
      --apply_images ${APPLY_LIST[*]}"
  else
    run "python \"$ANTSREG_PY\" \
      --fixed_image \"$ATLAS_FA\" \
      --moving_image \"$FA_REO\" \
      --output_prefix \"$OUT_PREFIX\" \
      --transform_type \"$TRANSFORM_TYPE\""
  fi

  # Rename applied outputs to space-Atlas (since we used apply_suffix "")
  # and keep the original metric naming; just add/replace space tag.
  # (Minimal & safe: only rename files that were just created from *_reoriented.)
  for src in "${MD_REO:-}" "${AD_REO:-}" "${RD_REO:-}"; do
    [[ -f "$src" ]] || continue
    b="$(basename "$src")"
    # strip _reoriented
    if [[ "$b" == *.nii.gz ]]; then
      base="${b%.nii.gz}"
      ext=".nii.gz"
    else
      base="${b%.nii}"
      ext=".nii"
    fi
    base="${base%_reoriented}"
    # enforce space-Atlas
    if [[ "$base" == *"_space-"* ]]; then
      outbase="$(echo "$base" | sed -E 's/_space-[^_]+/_space-Atlas/')"
    else
      outbase="${base}_space-Atlas"
    fi
    # applied file currently named after src (because suffix=""), so move it
    cur="$OUT_DWI_DIR/$b"
    [[ -f "$cur" ]] || continue
    run "mv -f \"$cur\" \"$OUT_DWI_DIR/${outbase}${ext}\""
  done

  log "Done sub-${SUB}/${sesdir}"
done

log "All done."
