# Pulseq SMS-EPI example: from acquisition to image reconstruction

[Overview](#overview)  
[Create the Pulseq sequence files](#create-the-pulseq-sequence-files)  
[Acquire fMRI data](#acquire-fmri-data)  
[Reconstruct images](#reconstruct-images)  
[Example data set](#example-data-set)  


## Overview

This folder contains MATLAB code that implements 
a complete vendor-agnostic workflow for acquiring and reconstructing
SMS-EPI data for functional MRI.

The workflow implementation provided here has been tested under Linux,
but should be portable to other platforms with minimal changes since it only
relies on MATLAB and Python code.


## Accessing an example data set

A resting-state fMRI data set in a volunteer is available upon request.


## Create the Pulseq sequence files

Follow the instructions in the 
[README.md](../sequence/README.md) file
in the [../sequence/](../sequence/) folder 
to create the following sequence files:

* **cal.seq**: EPI calibration scan
* **2d.seq**: 2D fully sampled slice-GRAPPA calibration scan
* **mb6.seq**: SMS-EPI fMRI scan (4 frames)

For GE users, the corresponding files **cal.tar**, **2d.tar**, and **mb6.tar** will also be created.

Use the default parameters, which will create an SMS-EPI scan with SMS factor 6,
2.4 mm isotropic resolution, TE/TR = 30/800 ms, and matrix size 90x90x60 (matching the ABCD protocol). 


## Acquire fMRI data

This sections contains information on how to run the .seq/.tar files created in the previous step on your scanner.

### Install the Pulseq interpreter

Interpreters for Siemens and GE are distributed freely within each 
respective research user community.
To get a copy of the interpreter code and usage instructions:  
**Siemens users**: Please contact Maxim Zaitsev (<maxim.zaitsev@uniklinik-freiburg.de>)  
**GE users**: Please contact Jon-Fredrik Nielsen (<jfnielse@umich.edu>).

General instructions for operating the interpreters can also be obtained from these contacts.

The interpreter just needs to be installed once -- once installed, it can be used to execute any Pulseq sequence.

### Prescribe the FOV

When running a Pulseq scan, the user is free to choose the center placement and rotation of the field of view (FOV);
most other sequence parameters such as field of view and matrix size are hardcoded into the .seq file 
and cannot be changed by the operator at prescription time.
This means that care must be taken to ensure that the slice-GRAPPA calibration scan (2d.seq) has the same FOV placement
as the fMRI scan (mb6.seq).
A simple way to do this is to **center the FOV at iso-center for all Pulseq scans**.

### Set the number of 'runs'

For the fMRI scan(s) (mb6.seq), set the `runs` parameter on the user interface to achieve the desired number of total frames.
Recall that mb6.seq contains 4 frames, and that the volume (frame) TR is 0.8 sec. 
So for example, for a 5-minute fMRI scan (376 frames), set `runs = 94`.

### Siemens-specific scan instructions

On Siemens, the .seq files can be executed in any order.

The data is written to disk in the usual .dat file format, as for any other sequence;
Pulseq does not alter this format in any way.
Therefore, the raw data file can be loaded with the usual Siemens data loading scripts.

### GE-specific scan instructions

On GE, the following scan receipe will ensure correct receive gain settings for all scans:

1. First load 'cal.seq', select auto prescan (APS), and run the scan.
2. Then load '2d.seq', enter manual prescan (MPS) and exit, and run the scan. 
3. Finally, load 'mb6.seq', enter MPS and exit, and run the scan.
4. Additional fMRI runs can now be executed by repeating the previous step,
   as long as the slice prescription on scanner (FOV placement) does not change.

The data is written to disk in the usual ScanArchive format, as for any other sequence.


## Reconstruct images

In the workflow demonstrated here, 
converting the raw (k-space) data to an image time-series involves the following steps:

1. Load the scanner-specific data file and convert to a custom HDF5 format. 
The reasons for doing this file format conversion are 
(1) the subsequent workflow is fully vendor-agnostic, and
(2) this custom HDF5 format splits the data frames across multiple smaller files
that can be loaded relatively quickly.
2. For each raw data set (cal.seq, 2d.seq, mb6.seq), 
   interpolate the ramp-sampled EPI data to a Cartesian grid using nufft.
3. Estimate odd/even EPI ghost correction parameters from the cal.seq scan, and apply this correction to the data.
4. Perform slice-GRAPPA calibration using data from the 2d.seq scan.
5. Reconstruct the data from the mb6.seq scan by performing slice-GRAPPA reconstruction for each slice group in each frame.

The script **recon_example.m** in the current folder implements these steps.
To use this script, do the following:
1. Edit **set_experimental_params.m**: Set the scanner vendor, data file paths, and other parameters.
2. Execute the reconstruction script:
   ```
   >> recon_example;
   ```

   See the comments in that script for further details.

