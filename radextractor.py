#!/usr/bin/env python3
import os
from pathlib import Path
import argparse
import nibabel as nib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from nilearn.image import resample_img

# ---------------------------
# Utility & QC Functions
# ---------------------------
def get_clean_name(path):
    return path.name.split('.nii')[0]

def generate_double_row_qc(data_img, roi_mask, wm_mask, output_path, title):
    """Generates a 6-subplot QC overlay: Top row = ROI, Bottom row = WM Mask."""
    data = data_img.get_fdata()
    roi_data = roi_mask.get_fdata()
    wm_data = wm_mask.get_fdata()
    
    coords = [s // 2 for s in data.shape]
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Planes for both rows
    planes = [
        (data[coords[0], :, :], "Sagittal"),
        (data[:, coords[1], :], "Coronal"),
        (data[:, :, coords[2]], "Axial")
    ]
    
    for i, (d_slice, name) in enumerate(planes):
        # Row 0: ROI Overlay (e.g., Corpus Callosum)
        axes[0, i].imshow(d_slice.T, cmap='gray', origin='lower')
        m_slice = [roi_data[coords[0], :, :], roi_data[:, coords[1], :], roi_data[:, :, coords[2]]][i]
        axes[0, i].imshow(m_slice.T, cmap='jet', alpha=0.3, origin='lower')
        axes[0, i].set_title(f"ROI: {name}")
        axes[0, i].axis('off')

        # Row 1: WM Mask Overlay
        axes[1, i].imshow(d_slice.T, cmap='gray', origin='lower')
        w_slice = [wm_data[coords[0], :, :], wm_data[:, coords[1], :], wm_data[:, :, coords[2]]][i]
        axes[1, i].imshow(w_slice.T, cmap='hot', alpha=0.3, origin='lower')
        axes[1, i].set_title(f"WM Mask: {name}")
        axes[1, i].axis('off')
        
    plt.suptitle(title, fontsize=16)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

# ---------------------------
# Core Logic
# ---------------------------
def main(args):
    bids_root = Path(args.bids_root)
    atlas_ver_dir = Path(args.atlas_dir) / f"{args.atlas_version}Atlas"
    
    # Paths
    roi_path = atlas_ver_dir / f"Regions_of_Interest_{args.atlas_version}" / "ROI_thr0"
    tpm_path = atlas_ver_dir / f"Tissue_Probability_Maps_{args.atlas_version}"
    
    # 1. Load Original Atlas Masks
    print(f"Loading atlas labels and TPMs...")
    roi_files = list(roi_path.glob("*.nii*"))
    raw_masks = {get_clean_name(f).replace("_thr0", ""): nib.load(str(f)) for f in roi_files}
    
    tpm_files = list(tpm_path.glob("*.nii*"))
    raw_tpms = {get_clean_name(f): nib.load(str(f)) for f in tpm_files}
    
    deriv_dir = bids_root / "derivatives" / "atlas_space-clean"
    subjects = [d for d in deriv_dir.iterdir() if d.is_dir() and d.name.startswith("sub-")]

    # Cache for resampled masks to avoid redundant work across sessions with the same grid
    current_resampled_rois = {}
    current_resampled_tpms = {}
    last_grid_params = None 

    for subj in subjects:
        for ses_dir in [s for s in subj.iterdir() if s.is_dir() and s.name.startswith("ses-")]:
            dwi_dir = ses_dir / "dwi"
            if not dwi_dir.exists(): continue
            
            fa_files = list(dwi_dir.glob("*desc-FA_dwi.nii*"))
            if not fa_files: continue
            
            ref_img = nib.load(str(fa_files[0]))
            grid_params = (ref_img.shape, np.sum(ref_img.affine))

            # 2. Resample ONLY if the grid has changed
            if grid_params != last_grid_params:
                print(f"Target grid change detected at {subj.name}. Resampling atlas...")
                current_resampled_rois = {}
                current_resampled_tpms = {}
                
                # Resample ROIs with requested flags
                for name, mask_img in raw_masks.items():
                    current_resampled_rois[name] = resample_img(
                        mask_img, target_affine=ref_img.affine, target_shape=ref_img.shape, 
                        interpolation='nearest', force_resample=True, copy_header=True
                    )
                # Resample TPMs with requested flags
                for name, tpm_img in raw_tpms.items():
                    current_resampled_tpms[name] = resample_img(
                        tpm_img, target_affine=ref_img.affine, target_shape=ref_img.shape, 
                        interpolation='nearest', force_resample=True, copy_header=True
                    )
                last_grid_params = grid_params

            csv_dir = ses_dir / "csv"
            csv_dir.mkdir(exist_ok=True, parents=True)

            # 3. 6-Subplot QC (CC ROI + WM TPM)
            qc_roi = current_resampled_rois.get("Corpus_Callosum")
            # Try to find a White Matter TPM (adjust string match if necessary)
            qc_wm = next((v for k, v in current_resampled_tpms.items() if "white" in k.lower() or "wm" in k.lower()), None)
            
            if qc_roi and qc_wm:
                qc_path = csv_dir / f"{subj.name}_{ses_dir.name}_Full_QC.png"
                generate_double_row_qc(ref_img, qc_roi, qc_wm, qc_path, f"QC Overlay: {subj.name} {ses_dir.name}")

            # 4. Extract Metrics
            for dti in ["FA", "MD", "AD", "RD"]:
                img_path = list(dwi_dir.glob(f"*desc-{dti}_dwi.nii*"))
                if not img_path: continue
                
                data_img = nib.load(str(img_path[0]))
                data = data_img.get_fdata()
                
                results = {}
                # Combine ROIs and TPMs for extraction
                combined_targets = {**current_resampled_rois, **current_resampled_tpms}
                
                for name, res_mask in combined_targets.items():
                    mask_data = res_mask.get_fdata()
                    # Apply thresholds and exclude background 0s
                    vals = data[(mask_data > 0.5) & (data > 0)]
                    
                    if vals.size == 0:
                        results[name] = [np.nan]*5
                    else:
                        results[name] = [np.nanmean(vals), np.nanstd(vals), np.nanmin(vals), np.nanmax(vals), np.median(vals)]
                
                df = pd.DataFrame(results, index=["mean", "std", "min", "max", "median"])
                df.index.name = "metric"
                df.to_csv(csv_dir / f"{dti}_metrics.csv")
                print(f"Saved metrics for {subj.name} {ses_dir.name} {dti}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bids-root", required=True)
    parser.add_argument("--atlas-dir", default="atlas")
    parser.add_argument("--atlas-version", default="4wk", choices=["4wk", "12wk"])
    main(parser.parse_args())
