# buildBIDS

Utilities for converting the internal HarmonizedMRI data organization into a BIDS-compatible directory structure for sharing with collaborators.

## Purpose

The scripts in this folder copy and rename imaging data from the internal archive into a standardized BIDS dataset. Where necessary, they also perform small modifications required for compatibility with downstream software.

Current functionality:

* Copy product BOLD NIfTI images into the BIDS `func/` directory.
* Copy Pulseq BOLD NIfTI images into the BIDS `func/` directory.
* Repair Pulseq NIfTI headers using the corresponding product image as a reference.
* Rename files using BIDS conventions.

Future versions will also support:

* MP-RAGE anatomical images (`anat/`)
* Dual-echo GRE field maps (`fmap/`)
* Raw ISMRMRD (`.mrd`) data (`sourcedata/`)
* BIDS sidecar JSON files
* Dataset-level files (`participants.tsv`, `dataset_description.json`, etc.)

## Expected source directory

Each scan session should be organized as

```text
srcRoot/
└── sub00012-umich-750MR-20250115-1/
    ├── product/
    │   ├── task_run1.h5.nii
    │   └── task_run2.h5.nii
    └── pulseq/
        ├── task_run1.h5.nii
        └── task_run2.h5.nii
```

The session directory name is parsed automatically to determine the BIDS subject and session labels.

## Output directory

The generated BIDS dataset is organized as

```text
bidsRoot/
└── sub-00012/
    └── ses-sub00012umich750MR202501151/
        └── func/
            ├── sub-00012_ses-sub00012umich750MR202501151_task-rest_acq-product_run-01_bold.nii
            └── sub-00012_ses-sub00012umich750MR202501151_task-rest_acq-pulseq_run-01_bold.nii
```

## Pulseq NIfTI processing

The reconstructed Pulseq images are modified before being written to the BIDS dataset.

The following operations are performed:

* Copy the NIfTI header from the corresponding product image.
* Preserve the reconstructed image dimensions.
* Set the voxel size to 2.4 × 2.4 × 2.4 mm.
* Set the repetition time (TR) to 0.8 s.
* Flip the first image dimension to match the product orientation.
* Scale the image to `int16`.
* Write the corrected NIfTI image.

## Usage

The current implementation operates on a single session.

```matlab
info = parseSessionName( ...
    'sub00012-umich-750MR-20250115-1', ...
    srcRoot, ...
    bidsRoot);

copyBOLD(info);
```

A future top-level `buildBIDS.m` script will iterate over all sessions and call the appropriate modality-specific functions.

## Notes

* Existing output files are not overwritten by default.
* Product images are copied without modification.
* Pulseq images use the corresponding product image as the reference for spatial orientation and NIfTI header information.
* The helper functions in the `private/` directory are intended to be shared by future conversion functions (e.g., `copyAnat`, `copyGRE`, etc.).

