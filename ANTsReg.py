#!/usr/bin/env python3
"""
ANTs Registration + Apply Transforms (lesion-aware, robust)
"""

import argparse
import os
import shutil
import ants
import sys

class AntsRegistration:
    def __init__(self, fixed_image_path, moving_image_path=None, output_prefix=None, moving_mask_path=None, transform_files=None):
        self.fixed_image_path = fixed_image_path
        self.moving_image_path = moving_image_path
        self.moving_mask_path = moving_mask_path
        self.output_prefix = output_prefix
        self.transform_files = transform_files
        self.transform = None

    def register(self, transform_type="SyN", reg_iterations=None):
        if self.moving_image_path is None and self.transform_files is None:
            raise ValueError("Must provide moving image or transform_files for registration.")

        fixed = ants.image_read(self.fixed_image_path)
        moving = ants.image_read(self.moving_image_path) if self.moving_image_path else None
        mask = ants.image_read(self.moving_mask_path) if self.moving_mask_path else None

        # Set default iterations
        if reg_iterations is None:
            if transform_type == "SyN":
                reg_iterations = (100, 50, 30)
            elif transform_type == "SyNAggro":
                reg_iterations = (100, 50, 30)

        # Run registration
        try:
            self.transform = ants.registration(
                fixed=fixed,
                moving=moving,
                type_of_transform=transform_type,
                moving_mask=mask,
                reg_iterations=reg_iterations,
                verbose=True
            )
        except Exception as e:
            raise RuntimeError(f"ANTs registration failed: {e}")

        # Check registration result
        if self.transform is None:
            raise RuntimeError(
                "ANTs registration returned None. "
                "Check: file paths, mask alignment, ANTsPy version, and extreme deformation."
            )

        # Save warped moving image
        if self.output_prefix:
            warped = self.transform.get("warpedmovout", None)
            if warped is not None:
                out_warped = f"{self.output_prefix}_warped.nii.gz"
                ants.image_write(warped, out_warped)
                print(f"Warped moving image saved: {out_warped}")

            # Save transforms
            fwd = self.transform.get("fwdtransforms", [])
            for i, tf in enumerate(fwd):
                shutil.copy(tf, f"{self.output_prefix}_fwd_{i}_{os.path.basename(tf)}")

            inv = self.transform.get("invtransforms", [])
            for i, tf in enumerate(inv):
                shutil.copy(tf, f"{self.output_prefix}_inv_{i}_{os.path.basename(tf)}")

    def apply_transform_one(self, image_path, output_path, is_label=False):
        if self.transform:
            transformlist = self.transform["fwdtransforms"]
        elif self.transform_files:
            transformlist = self.transform_files
        else:
            raise ValueError("Run registration first or provide transform_files.")

        fixed = ants.image_read(self.fixed_image_path)
        moving = ants.image_read(image_path)
        interp = "nearestNeighbor" if is_label else "linear"

        out = ants.apply_transforms(
            fixed=fixed,
            moving=moving,
            transformlist=transformlist,
            interpolator=interp
        )
        ants.image_write(out, output_path)
        print(f"Transformed image saved: {output_path}")

    def apply_transform_many(self, image_paths, output_dir, suffix="_space-Atlas", is_label=False):
        os.makedirs(output_dir, exist_ok=True)
        for p in image_paths:
            base = os.path.basename(p)
            if base.endswith(".nii.gz"):
                out_base = base[:-7] + suffix + ".nii.gz"
            elif base.endswith(".nii"):
                out_base = base[:-4] + suffix + ".nii"
            else:
                raise ValueError(f"Unsupported file extension: {p}")
            out_path = os.path.join(output_dir, out_base)
            self.apply_transform_one(p, out_path, is_label=is_label)


def main():
    parser = argparse.ArgumentParser(description="ANTsPy Registration + Apply Transforms")
    parser.add_argument("--fixed_image", type=str, required=True)
    parser.add_argument("--moving_image", type=str)
    parser.add_argument("--moving_mask", type=str, help="Lesion mask (binary) for moving image")
    parser.add_argument("--output_prefix", type=str, required=True)
    parser.add_argument("--transform_type", type=str, default="SyN")
    parser.add_argument("--reg_iterations", nargs=3, type=int, help="Iterations for registration")

    parser.add_argument("--transform_files", nargs="+", help="Precomputed transforms")
    parser.add_argument("--apply_images", nargs="*", default=[], help="Images to apply transform to")
    parser.add_argument("--apply_out_dir", type=str, help="Output dir for applied images")
    parser.add_argument("--apply_suffix", type=str, default="_space-Atlas")
    parser.add_argument("--is_label", action="store_true", help="Nearest-neighbor interpolation for labels")

    args = parser.parse_args()

    reg = AntsRegistration(
        fixed_image_path=args.fixed_image,
        moving_image_path=args.moving_image,
        moving_mask_path=args.moving_mask,
        output_prefix=args.output_prefix,
        transform_files=args.transform_files
    )

    # Run registration
    if args.moving_image:
        reg.register(transform_type=args.transform_type, reg_iterations=tuple(args.reg_iterations) if args.reg_iterations else None)

    # Apply to images
    if args.apply_images:
        if not args.apply_out_dir:
            sys.exit("--apply_out_dir is required when using --apply_images")
        reg.apply_transform_many(
            image_paths=args.apply_images,
            output_dir=args.apply_out_dir,
            suffix=args.apply_suffix,
            is_label=args.is_label
        )


if __name__ == "__main__":
    main()
