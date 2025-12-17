#!/usr/bin/env python3
import argparse
import os
import numpy as np
import nibabel as nib


def parse_int_list(s: str):
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def sanitize_array(a: np.ndarray) -> np.ndarray:
    a = np.nan_to_num(a, nan=0.0, posinf=0.0, neginf=0.0)
    return np.ascontiguousarray(a)


def squeeze_singleton_4th(a: np.ndarray) -> np.ndarray:
    if a.ndim == 4 and a.shape[-1] == 1:
        return np.ascontiguousarray(a[..., 0])
    return a


def out_with_suffix(in_path: str, suffix: str) -> str:
    if in_path.endswith(".nii.gz"):
        return in_path[:-7] + suffix + ".nii.gz"
    if in_path.endswith(".nii"):
        return in_path[:-4] + suffix + ".nii"
    raise SystemExit(f"Unsupported extension (need .nii or .nii.gz): {in_path}")


def reorient_voxel_only(img: nib.Nifti1Image, permute_axes, flip_axes):
    data = img.get_fdata(dtype=np.float32)
    data = sanitize_array(data)
    data = squeeze_singleton_4th(data)

    if data.ndim != 3:
        raise SystemExit(f"Expected 3D after squeeze, got shape {data.shape}")

    if permute_axes is not None:
        if sorted(permute_axes) != [0, 1, 2]:
            raise SystemExit(f"--permute_axes must be permutation of 0,1,2. Got {permute_axes}")
        data = np.ascontiguousarray(np.transpose(data, tuple(permute_axes)))

    if flip_axes:
        for ax in flip_axes:
            if ax not in (0, 1, 2):
                raise SystemExit(f"--flip_axes values must be 0/1/2. Got {ax}")
            data = np.ascontiguousarray(np.flip(data, axis=ax))

    # IMPORTANT: keep affine unchanged (no compensation) so physical orientation changes
    out = nib.Nifti1Image(data.astype(np.float32), img.affine)
    out.header.set_data_dtype(np.float32)
    return out


def main():
    p = argparse.ArgumentParser(
        description="Reorient NIfTI by voxel permute+flip (nibabel). Keeps affine unchanged (NO compensation)."
    )
    p.add_argument("--in", dest="in_path", required=True, help="Input NIfTI")
    p.add_argument("--out", dest="out_path", default=None, help="Output NIfTI. If omitted, writes *_reoriented.")
    p.add_argument("--permute_axes", default=None, help='Permutation like "0,2,1"')
    p.add_argument("--flip_axes", default=None, help='Flip axes like "2" or "0,2"')
    p.add_argument("--canonicalize", action="store_true", help="nib.as_closest_canonical before voxel ops")
    p.add_argument("--suffix", default="_reoriented", help="Suffix if --out not provided (default: _reoriented)")
    args = p.parse_args()

    perm = parse_int_list(args.permute_axes) if args.permute_axes else None
    flips = parse_int_list(args.flip_axes) if args.flip_axes else []

    img = nib.load(args.in_path)
    if args.canonicalize:
        img = nib.as_closest_canonical(img)

    out_img = reorient_voxel_only(img, perm, flips)

    out_path = args.out_path or out_with_suffix(args.in_path, args.suffix)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    nib.save(out_img, out_path)


if __name__ == "__main__":
    main()
