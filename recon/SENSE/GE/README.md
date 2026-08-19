# Preparing SMS-EPI Data for Reconstruction (GE Scanner)

This README explains how to prepare for reconstruction
data acquired on a GE scanner
using a SMS-EPI scan
created with TOPPE.

## Dependencies

First,
make sure you have downloaded
this [SMS-EPI](https://github.com/HarmonizedMRI/SMS-EPI)
repository with, e.g., `git clone`:
```bash
$ git clone https://github.com/HarmonizedMRI/SMS-EPI
```

This data preparation pipeline requires
a working MATLAB installation.
The instructions in this README
were tested using R2022b Update 3,
though other versions likely would still work.

The following MATLAB toolboxes are needed:

1. [`toppe`](https://github.com/toppeMRI/toppe),
1. [`mirt`](https://github.com/JeffFessler/mirt),
1. [`hmriutils`](https://github.com/HarmonizedMRI/utils), and
1. [`bart`](https://mrirecon.github.io/bart/installation.html).

For `toppe` and `hmriutils`,
you can download them
with, e.g., `git clone`:
```bash
$ git clone https://github.com/toppeMRI/toppe
$ git clone https://github.com/JeffFessler/mirt
$ git clone https://github.com/HarmonizedMRI/utils
```
For `bart`,
follow the installation instructions
provided by the above link.

After following the above steps,
everything should be set up
for preparing SMS-EPI data for reconstruction.

## Step 0: Preliminaries

Before starting,
navigate to
the directory
that will serve as your workspace
for the following.

First,
you will need the GE data files (P-files)
that contain the scan data that was acquired.
Typically,
there will be two such files:
one containing data
that will be used for computing sensitivity maps,
and the other containing the SMS-EPI data.
In this README,
the two P-files are
`P_b0.7`
and `P_fmri_mb2.7`,
respectively.
(`P_b0.7` is so named
because it contains two spoiled gradient recalled echo (SPGR) scans
that can be used for B0 field mapping.
The first of those scans, however,
can be used for computing sensitivity maps.
And `P_fmri_mb2.7` is so named
because the data was acquired
using a fMRI protocol
with a multiband (mb) factor of 2.)
P-files can be downloaded
from [Deep Blue Data](TODO)
(link pending)
for those interested
in following along
with this README.

Next,
you will need the TOPPE scan files.
These are typically grouped into .tar files,
one per scan.
In this README,
the two scan files are
`b0.tar`
and `fmri_mb2.tar`.
Actually,
you really only need one file
from each .tar file.
From `b0.tar` you just need `readout.mod`,
and from `fmri_mb2.tar` you just need `epiInfo.mat`.
You can extract them with, e.g.,
```bash
$ tar xf b0.tar readout.mod
$ tar xf fmri_mb2.tar epiInfo.mat
```
from the command line.
The .tar files
associated with the previously mentioned P-files
can also be downloaded
from [Deep Blue Data](TODO).

Finally,
create a MATLAB script
that sets up the depenencies
for your current session.
The MATLAB script should contain
```matlab
% Include the toppe namespace
addpath path/to/toppe/ % TODO: Change to the correct path

% Set up MIRT
cwd = cd('path/to/mirt/'); % TODO: Change to the correct path
setup;
cd(cwd);

% Include the hmriutils namespace
addpath path/to/HarmonizedMRI/utils/ % TODO: Change to the correct path

% Set up bart
setenv('TOOLBOX_PATH', 'path/to/bart'); % TODO: Change to the correct path
addpath(strcat(getenv('TOOLBOX_PATH'), '/matlab'));

% Provide access to functions in the SMS-EPI repository
addpath path/to/HarmonizedMRI/SMS-EPI/recon/GE/ % TODO: Change to the correct path
```
(make sure to update the paths).
In this README,
this script is called `setup.m`.

After collecting all the necessary files,
your working directory should contain
(at a minimum)
the following files:
```
P_b0.7
P_fmri_mb2.7
readout.mod
epiInfo.mat
setup.m
```

## Step 1: Loading and Reshaping the SMS-EPI Data

Start a MATLAB session,
`cd` to your working directory,
and then run
```matlab
>> setup;
```
Next,
call `readframe`,
passing in the appropriate inputs:
```matlab
>> frame = 10; % Process frame 10 of the data (change this to whatever you want)
>> mb = 2; % In this README, the SMS multiband factor is 2
>> readframe(frame, mb, 'pfile', 'P_fmri_mb2.7', 'epiInfoFile', 'epiInfo.mat');
```
The call to `readframe` creates
a file called `frame010.mat`
containing the k-space data
and CAIPI sampling mask
for frame 10 of the SMS-EPI data.
`readframe` also reshapes the data
into an appropriate format
for reconstruction.

## Step 2: Computing the Sensitivity Maps

In the same MATLAB session,
call `getsense`,
passing in the appropriate inputs:
```matlab
>> getsense('pfile', 'P_b0.7', 'readoutFile', 'readout.mod', 'slices', 7:66);
```
(The `'slices'` keyword argument
is specified because `P_b0.7` contains data
for 72 z phase encodes
to account for slab profile effects,
but the SMS-EPI data
contains only 60 slices.
Thus we want to compute sensitivity maps
only for the middle 60 slices.)
The call to `getsense` creates
a file called `sense.mat`
containing the sensitivity maps.
Note that `getsense` can take a while to finish.

## Next Steps

Now that the SMS-EPI data
has been loaded and reshaped
and the sensitivity maps
have been estimated,
you can proceed to reconstruct the images
using the [SENSE-based reconstruction method in this repository](..).
