# Reconstructing SMS-EPI Data

## SENSE Reconstruction

[`sense.jl`](sense.jl) defines the function `recon`
for reconstructing a 3D volume from SMS-EPI data
using a SENSE-based image reconstruction method.
This section of the README
walks through how to use `recon`.

### Dependencies

First,
make sure you have downloaded
this [SMS-EPI](https://github.com/HarmonizedMRI/SMS-EPI)
repository with, e.g., `git clone`:
```bash
$ git clone https://github.com/HarmonizedMRI/SMS-EPI
```

This reconstruction pipeline requires
a working Julia installation.
The latest version is recommended,
though earlier versions likely would still work.
You can [download the latest Julia version here](https://julialang.org/downloads/#current_stable_release).

`recon` depends on two Julia packages:
[MAT.jl](https://github.com/JuliaIO/MAT.jl) and
[EPI3D.jl](https://github.com/HarmonizedMRI/3DEPI/tree/main/recon/sense),
and `showrecon` depends on one additional Julia package:
[MIRTjim.jl](https://github.com/JeffFessler/MIRTjim.jl).
In order to use EPI3D.jl,
you must first download or clone the repository, e.g.,
```bash
$ git clone https://github.com/HarmonizedMRI/3DEPI
```
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
(@v1.9) pkg> activate .

(working_dir) pkg> add MAT MIRTjim

(working_dir) pkg> dev path/to/3DEPI/recon/sense
```
(Make sure to change `path/to/`
to the location where you downloaded or cloned
the 3DEPI repository.)

### Step 0: Preliminaries

Before starting,
navigate to
the directory
that will serve as your workspace
for the following.
(This directory should be the same
as the one used above
to set up the Julia package environment.)

First,
you will need files
that contain the input
`recon` expects.
Specifically,
`recon` needs the SMS-EPI k-space data `data`,
a mask `mask` indicating the k-space locations
of each of the data samples
acquired using a CAIPI sampling scheme,
and the coil sensitivity maps `smap`.
These inputs can be given directly to `recon`,
or they can be stored in .mat files
to be loaded in `recon`.
In the latter case,
two .mat files are needed:
one that stores `data` and `mask`,
and another that stores `smap`.
In this README,
the two .mat files are
`frame010.mat`
and `sense.mat`.
(These files can be generated, e.g.,
by following the directions
in the [README for preparing GE SMS-EPI data for reconstruction](GE/README.md),
or can be downloaded
from [Deep Blue Data](TODO).)

Finally,
create a Julia script
that sets up the depenencies
for your current session.
The Julia script should contain
```julia
# Make sure the package environment is active
import Pkg
Pkg.activate(@__DIR__)

# Load functions in the SMS-EPI repository
include("path/to/HarmonizedMRI/SMS-EPI/recon/sense.jl") # TODO: Change to the correct path
include("path/to/HarmonizedMRI/SMS-EPI/recon/showrecon.jl") # TODO: Change to the correct path
```
(make sure to update the paths).
In this README,
this script is called `setup.jl`.

After collecting all the necessary files,
your working directory should contain
(at a minimum)
the following files:
```
frame010.mat
sense.mat
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

### Step 1: Reconstruction

Start a Julia session
in your working directory
(or `cd` to your working directory
after starting Julia),
and then run
```julia
include("setup.jl")
```
(Note that the above `include`
needs to be done just once per Julia session.)
Next,
call `recon`,
passing in the appropriate inputs:
```julia
recon("frame010.mat", "sense.mat"; lambda = 50, outfile = "recon010.mat")
```
The call to `recon` creates
a file called `recon010.mat`
containing the reconstructed image.
See the documentation for a description of the inputs to `recon`.

### Step 2: Visualization of the Reconstructed Image

The following will display the reconstructed image:
```julia
showrecon("recon010.mat"; title = "SMS-EPI Reconstruction (mb = 2)", yflip = false)
```
See the documentation for a description of the inputs to `showrecon`.

### Documentation

See the docstrings in [`sense.jl`](sense.jl)
and [`showrecon.jl`](showrecon.jl),
or run the following from the Julia REPL
after running `include("setup.jl")`:
```julia-repl
julia> ? # enter Julia's help mode
help?> recon

help?> showrecon
```
