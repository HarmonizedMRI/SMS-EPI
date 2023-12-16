# A complete SMS-EPI example: from acquisition to image reconstruction

This folder contains MATLAB code that implements 
a complete **vendor-agnostic workflow for acquiring and reconstructing
SMS-EPI data for functional MRI**.
This is part of the HarmonizedMRI project for reproducible MRI, 
see https://harmonizedmri.github.io/.

The SMS-EPI fMRI pulse sequence is defined in the 
[Pulseq](https://pulseq.github.io/ "Vendor-agnostic MRI pulse sequences")
vendor-agnostic file format,
that allows the same sequence to be executed on multiple MRI
hardware platforms; at present, all Siemens and GE scanner models 
are supported.
The fMRI image time-series is reconstructed using slice GRAPPA.
In this way, you have complete control over, and knowledge of, 
the entire protocol from acquisition to reconstructed images;
no vendor-dependent reconstruction or image processing is performed.

The SMS-EPI fMRI sequence provided here matches the ABCD protocol
(2.4 mm isotropic resolution, 90x90x60 matrix size, 0.8s temporal resolution,
SMS factor 6), and has been tested extensively in phantoms and volunteers;
however you are free the change these parameters as described below.
You may also want to add your own image post-processing steps;
here we simply guarantee that the SMS-EPI images are acquired and reconstructed
in a fully transparent, vendor-agnostic, and customizable way.

The workflow involves the following steps:

## Step 1: Installation

The code in this repository can be
downloaded using the following `git clone` command 
from, e.g., a Linux shell or Git Bash:
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

**Linux users** can run `set_path.m` 
to perform these actions entirely from within the MATLAB prompt, e.g.:
```
>> system('git clone --branch develop git@github.com:HarmonizedMRI/SMS-EPI.git');
>> cd SMS-EPI/example
>> set_paths;
```

**Windows users**: for instructions on how to do `git clone` in Windows, 
see [this video](https://www.youtube.com/watch?v=Av7lcVIbEBY&t=1s).

**MAC users?**

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
in MATLAB do:
```
>> write_sequences;
```
This creates the following .seq files:
* `sms-epi.seq`: SMS-EPI scan for fMRI
* `cal.seq`: EPI ghost calibration scan
* `acs.seq`: individual slice acquisitions for slice GRAPPA calibration

TODO: blip up/down spin-echo scan for distortion correction

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

## TODO  
   * blip up/down spin-echo scan for distortion correction
   * Siemens example data
   * In vivo example data (3D GRE ACS scan)
   * Share data on DeepBlue/other
   * Consider using similarity b/w calibration scans and recon'd images as QC metric


