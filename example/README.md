# A complete Pulseq SMS-EPI example: from acquisition to image reconstruction

This folder contains MATLAB code that implements 
a complete vendor-agnostic workflow for acquiring and reconstructing
SMS-EPI data for functional MRI.

The specific workflow implementation provided here is designed for **Linux**;
If you'd like to implement a similar workflow on other operating systems,
let us know and we'll work with you to make that happen.


## Step 1: Create the Pulseq sequence files

Follow the instructions in the 
[README.md](../sequence/README.md) file
in the [../sequence/](../sequence/) folder 
to create the following sequence files:

* **cal.seq**: EPI calibration scan
* **2d.seq**: 2D fully sampled slice-GRAPPA calibration scan
* **mb6.seq**: SMS-EPI fMRI scan (4 frames)

For GE users, the corresponding files **cal.tar**, **2d.tar**, and **mb6.tar** will also be created.


## Step 2: Acquire fMRI data

This sections contains information on how to run the .seq/.tar files created in Step 1 on your scanner.

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

For the fMRI scan(s) (mb6.seq), set the 'runs' parameter on the user interface to achieve the desired number of total frames.
Recall that mb6.seq contains 4 frames, and that the volume (frame) TR is 0.8 sec. 
So for example, for a 5-minute fMRI scan (376 frames), set 'runs' to 94.

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


## Reconstruct the fMRI image time-series

## Step 5: Reconstruct time-series images
Place the acquired data files in the folder specified
in `set_experimental_params.m`. 
Then execute the following MATLAB commands:
```
set_experimental_params;       % get experimental parameters
get_ghost_calibration_data;    % get EPI ghost calibration data 
get_acs_data;                  % get individual slice (ACS) data for slice GRAPPA
recon_timeseries;              % do slice GRAPPA reconstruction
```



# OLD

This folder contains MATLAB code that implements 
a complete vendor-agnostic workflow for acquiring and reconstructing
SMS-EPI data for functional MRI.

The specific workflow implementation provided here is 
designed for **Linux**;
If you'd like to implement a similar workflow on
other operating systems, let us know and we'll work with
you to make that happen.

The SMS-EPI fMRI pulse sequence is defined in the 
[Pulseq](https://pulseq.github.io/ "Vendor-agnostic MRI pulse sequences")
vendor-agnostic file format,
that allows the same sequence definition to be executed on multiple MRI
hardware platforms; at present, Siemens and GE scanner models 
are supported.
The fMRI image time-series is reconstructed using slice GRAPPA.
In this way, you have complete control over, and knowledge of, 
the entire protocol from acquisition to reconstructed images;
no vendor-dependent reconstruction or image processing is performed.

The SMS-EPI fMRI sequence provided here matches the ABCD protocol
(2.4 mm isotropic resolution, 90x90x60 matrix size, 0.8s temporal resolution,
SMS factor 6), and has been tested in phantoms and volunteers;
however you are free to change these parameters as described below.
You may also want to add your own image post-processing steps;
here we simply guarantee that the SMS-EPI images are acquired and reconstructed
in a fully transparent, vendor-agnostic, and customizable way.

The workflow involves the following steps:

## Step 0: Requirements

**MATLAB** is required to run the code in this repository.

**git** is the most convenient way to download the code as described in the next step. 
Alternatively, you can navigate to each repository and download the code as a ZIP file 
(under the green "Code" button in each repository).

## Step 1: Installation

The code in this repository can be
downloaded using the following `git clone` command 
from, e.g., a Linux shell or a Git Bash shell in Windows:
```
git clone --branch develop git@github.com:HarmonizedMRI/SMS-EPI.git
```

In addition, the code dependencies must be downloaded:
```
git clone git@github.com:pulseq/pulseq.git
git clone --branch develop git@github.com:HarmonizedMRI/utils.git
git clone --branch develop git@github.com:HarmonizedMRI/PulCeq.git
git clone --branch develop git@github.com:toppeMRI/toppe.git
git clone git@github.com:JeffFessler/mirt.git
```
followed by the corresponding `addpath` commands in MATLAB:
```
addpath pulseq/matlab      % Pulseq (+mr namespace)
addpath utils              % +hmriutils namespace 
addpath PulCeq/matlab      % only needed if converting to GE scan files
addpath toppe              % only needed if converting to GE scan files
cd mirt; setup; cd ..;     % 'im' function for display
```

**Linux users** can run `get_code_and_set_paths.m` 
to perform these actions entirely from within the MATLAB prompt, e.g.:
```
>> system('git clone --branch develop git@github.com:HarmonizedMRI/SMS-EPI.git');
>> cd SMS-EPI/example
>> get_code_and_set_paths;
```

**Windows users**: for instructions on how to do `git clone` in Windows, 
see [this video](https://www.youtube.com/watch?v=Av7lcVIbEBY&t=1s).

## Step 2: Set acquisition and other experimental parameters
You may edit the file `set_experimental_params.m` to set things like
spatial resolution and the name and  location of your acquired data files.
However, when first testing this code, we recommend that you leave the acquisition
parameters as they are. Then execute it:
```
>> set_experimental_params;
```

## Step 3: Create the Pulseq (.seq) files

To create the SMS-EPI .seq file as well as the two calibration sequences, 
or to obtain a .seq file already created by the HarmonizedMRI study team,
contact jfnielse@umich.edu along with your Github user name
(this helps us track usage).

## Step 4: Run the Pulseq files on your scanner
This requires that you install the Pulseq interpreter for your scanner make and model.
Interpreters for Siemens and GE are distributed freely within each 
respective research user community.
To get a copy of the interpreter code and usage instructions:  
**Siemens users**: Please contact Maxim Zaitsev (<maxim.zaitsev@uniklinik-freiburg.de>)  
**GE users**: Please contact Jon-Fredrik Nielsen (<jfnielse@umich.edu>)  

## Step 5: Reconstruct time-series images
Place the acquired data files in the folder specified
in `set_experimental_params.m`. 
Then execute the following MATLAB commands:
```
set_experimental_params;       % get experimental parameters
get_ghost_calibration_data;    % get EPI ghost calibration data 
get_acs_data;                  % get individual slice (ACS) data for slice GRAPPA
recon_timeseries;              % do slice GRAPPA reconstruction
```

## Example data files 

Example phantom data will be posted here soon.


