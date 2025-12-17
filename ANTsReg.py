#!/usr/bin/env python3
import argparse
import os
import shutil
import ants


class AntsRegistration:
    def __init__(self, fixed_image_path, moving_image_path=None, output_prefix=None, transform_files=None):
        self.fixed_image_path = fixed_image_path
        self.moving_image_path = moving_image_path
        self.output_prefix = output_prefix
        self.transform_files = transform_files
        self.transform = None

    def register(self, transform_type="SyN", clean=False):
        if self.moving_image_path is None:
            raise ValueError("Moving image path must be provided for registration.")

        fixed_image = ants.image_read(self.fixed_image_path)
        moving_image = ants.image_read(self.moving_image_path)

        params = {}
        if transform_type == "SyN":
            params['reg_iterations'] = (100, 50, 30)

        self.transform = ants.registration(
            type_of_transform=transform_type,
            fixed=fixed_image,
            moving=moving_image,
            verbose=True,
            **params
        )

        if clean or not self.output_prefix:
            return

        warped = self.transform.get("warpedmovout", None)
        if warped is not None:
            ants.image_write(warped, f"{self.output_prefix}_warped.nii.gz")

        # Save forward transforms (useful for downstream)
        fwd = self.transform.get("fwdtransforms", None)
        if isinstance(fwd, list):
            for i, tf in enumerate(fwd):
                shutil.copy(tf, f"{self.output_prefix}_fwd_{i}_{os.path.basename(tf)}")

        # Bugfix: correct list check for invtransforms
        inv = self.transform.get("invtransforms", None)
        if isinstance(inv, list):
            for i, tf in enumerate(inv):
                shutil.copy(tf, f"{self.output_prefix}_inv_{i}_{os.path.basename(tf)}")

    def apply_transform_one(self, image_path, output_path, is_label=False):
        if self.transform:
            transformlist = self.transform["fwdtransforms"]
        elif self.transform_files:
            transformlist = self.transform_files
        else:
            raise ValueError("Run registration first or provide --transform_files.")

        fixed_image = ants.image_read(self.fixed_image_path)
        image = ants.image_read(image_path)

        interpolation = "nearestNeighbor" if is_label else "linear"
        out = ants.apply_transforms(
            fixed=fixed_image,
            moving=image,
            transformlist=transformlist,
            interpolator=interpolation
        )
        ants.image_write(out, output_path)

    def apply_transform_many(self, image_paths, output_dir, suffix="_space-Atlas", is_label=False):
        os.makedirs(output_dir, exist_ok=True)
        for p in image_paths:
            base = os.path.basename(p)
            if base.endswith(".nii.gz"):
                out_base = base[:-7] + suffix + ".nii.gz"
            elif base.endswith(".nii"):
                out_base = base[:-4] + suffix + ".nii"
            else:
                raise ValueError(f"Unsupported extension for {p}")

            out_path = os.path.join(output_dir, out_base)
            self.apply_transform_one(p, out_path, is_label=is_label)


def main():
    parser = argparse.ArgumentParser(description="ANTsPy registration + apply transforms (no reorientation).")
    parser.add_argument("--fixed_image", type=str, required=True, help="Fixed image")
    parser.add_argument("--moving_image", type=str, help="Moving image")
    parser.add_argument("--output_prefix", type=str, help="Prefix for registration outputs")
    parser.add_argument("--transform_type", type=str, default="SyN", help="ANTs transform type")

    parser.add_argument("--transform_files", nargs="+",
                        help="Transform files for apply (order: warps then affine). If provided, skips needing registration.")

    parser.add_argument("--apply_images", nargs="*", default=[],
                        help="One or more images to apply the (just-computed) transform to.")
    parser.add_argument("--apply_out_dir", type=str, default=None,
                        help="Output directory for applied images (required if --apply_images is used).")
    parser.add_argument("--apply_suffix", type=str, default="_space-Atlas",
                        help="Suffix appended to applied outputs (default: _space-Atlas)")
    parser.add_argument("--is_label", action="store_true", help="Use nearest-neighbor interpolation")

    args = parser.parse_args()

    reg = AntsRegistration(
        fixed_image_path=args.fixed_image,
        moving_image_path=args.moving_image,
        output_prefix=args.output_prefix,
        transform_files=args.transform_files
    )

    # If moving_image is provided, run registration
    if args.moving_image:
        reg.register(transform_type=args.transform_type, clean=(args.output_prefix is None))

    # Apply to many right away
    if args.apply_images:
        if not args.apply_out_dir:
            raise SystemExit("--apply_out_dir is required when using --apply_images")
        reg.apply_transform_many(
            image_paths=args.apply_images,
            output_dir=args.apply_out_dir,
            suffix=args.apply_suffix,
            is_label=args.is_label
        )


if __name__ == "__main__":
    main()
