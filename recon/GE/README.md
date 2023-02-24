# Reconstructing SMS-EPI Data (GE Scanner)

This README explains how to reconstruct data
acquired on a GE scanner
using a SMS-EPI scan
created with TOPPE.

## Dependencies

First,
make sure you have downloaded
this [SMS-EPI](https://github.com/HarmonizedMRI/SMS-EPI)
repository with, e.g., `git clone`.

This reconstruction pipeline requires
working MATLAB
and Julia installations.
The latest versions are recommended,
though earlier versions likely would still work.
You can download the latest Julia version
[here](https://julialang.org/downloads/#current_stable_release).

The following toolboxes are needed
for the MATLAB portion
of the processing pipeline:

1. [`toppe`](https://github.com/toppeMRI/toppe),
1. [`hmriutils`](https://github.com/HarmonizedMRI/utils), and
1. [`bart`](https://mrirecon.github.io/bart/installation.html).

For `toppe` and `hmriutils`,
you can download them
with, e.g., `git clone`.
For `bart`,
follow the installation instructions
provided by the above link.

Use `git clone`
to download the following code,
needed for the Julia portion
of the processing pipeline:

1. [EPI3D.jl](https://github.com/HarmonizedMRI/3DEPI).

To set up the Julia package environment,
first start Julia
in the directory
in which you will be working.
(If needed,
create such a directory first.)
Then set up dependencies
by running the following:
```julia-repl
julia> ] # enter Julia's package prompt
(@v1.8) pkg> activate .

(working_dir) pkg> add MAT MIRTjim

(working_dir) pkg> dev path/to/3DEPI/recon/sense
```
(Make sure to change `path/to/`
to the location where you cloned
the 3DEPI repository.)

After following the above steps,
everything should be set up
for processing SMS-EPI data.

## Step 0: Preliminaries

Before starting,
navigate to
the directory
that will serve as your workspace
for the following.
(This directory should be the same
as the one used above
to set up the Julia package environment.)

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
`P,b0.7`
and `P,fmri_mb6.7`,
respectively.
(`P,b0.7` is so named
because it contains two spoiled gradient recalled echo (SPGR) scans
that can be used for B0 field mapping.
The first of those scans, however,
can be used for computing sensitivity maps.
And `P,fmri_mb6.7` is so named
because the data was acquired
using a fMRI protocol
with a multiband (mb) factor of 6.)

Next,
you will need the TOPPE scan files.
These are typically grouped into .tar files,
one per scan.
In this README,
the two scan files are
`b0.tar`
and `fmri_mb6.tar`.
Actually,
you really only need one file
from each .tar file.
From `b0.tar` you just need `readout.mod`,
and from `fmri_mb6.tar` you just need `epiInfo.mat`.
You can extract them with, e.g.,
```
$ tar xf b0.tar readout.mod
$ tar xf fmri_mb6.tar epiInfo.mat
```
from the command line.

Finally,
create a MATLAB script
and a Julia script
that set up the depenencies
for your current session.
The MATLAB script should contain
```matlab
% Include the toppe namespace
addpath path/to/toppe/ % TODO: Change to the correct path

% Include the hmriutils namespace
addpath path/to/HarmonizedMRI/utils/ % TODO: Change to the correct path

% Set up bart
setenv('TOOLBOX_PATH', 'path/to/bart'); % TODO: Change to the correct path
addpath(strcat(getenv('TOOLBOX_PATH'), '/matlab'));

% Provide access to functions in the SMS-EPI repository
addpath path/to/HarmonizedMRI/SMS-EPI/recon/GE/ % TODO: Change to the correct path
```
(make sure to update the paths),
and the Julia script should contain
```julia
# Make sure the package environment is active
import Pkg
Pkg.activate(@__DIR__)

# Load functions in the SMS-EPI repository
include("path/to/HarmonizedMRI/SMS-EPI/recon/sense.jl") # TODO: Change to the correct path
include("path/to/HarmonizedMRI/SMS-EPI/recon/showrecon.jl") # TODO: Change to the correct path
```
(make sure to update the path).
In this README,
these two scripts are called
setup.m and setup.jl, respectively.

After collecting all the necessary files,
your working directory should contain
(at a minimum)
the following files:
```
P,b0.7
P,fmri_mb6.7
readout.mod
epiInfo.mat
setup.m
setup.jl
Project.toml
Manifest.toml
```
(Note that the .toml files
were created by the Julia package manager
when setting up the dependencies.
If these files are missing,
go to the Dependencies section
of this README
and follow the instructions
on how to set up the Julia package environment.)

## Step 1: Loading and Reshaping the Data (MATLAB)

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
>> mb = 6; % In this README, the SMS multiband factor is 6
>> readframe(frame, mb, 'pfile', 'P,fmri_mb6.7', 'epiInfoFile', 'epiInfo.mat');
```
The call to `readframe` creates
a file called `frame010.mat`
containing the k-space data
and CAIPI sampling mask
for frame 10 of the SMS-EPI data.
`readframe` also reshapes the data
into an appropriate format
for reconstruction.

## Step 2: Computing the Sensitivity Maps (MATLAB)

In the same MATLAB session,
call `getsense`,
passing in the appropriate inputs:
```matlab
>> getsense('pfile', 'P,b0.7', 'readoutFile', 'readout.mod');
```
The call to `getsense` creates
a file called `sense.mat`
containing the sensitivity maps.
Note that `getsense` can take a while to finish.

## Step 3: Reconstructing SMS-EPI Data (Julia)

Start a Julia session
in your working directory
(or `cd` to your working directory
after starting Julia),
and then run
```julia-repl
julia> include("setup.jl")
```
Next,
call `recon`,
passing in the appropriate inputs:
```julia-repl
julia> recon("frame010.mat", "sense.mat"; outfile = "recon010.mat");
```
The call to `recon` creates
a file called `recon010.mat`
containing the reconstructed image.

## Step 4: Visualizing the Reconstructed Data (Julia)

Now that the data has been reconstructed,
you can plot it.
You can do this
in whatever environment you desire,
but this README will cover
how to do so in Julia.
It is as simple as calling `showrecon`:
```julia-repl
julia> showrecon("recon010.mat"; title = "SMS-EPI Reconstruction (mb = 6)")
```
