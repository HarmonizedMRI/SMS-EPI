# Reconstructing SMS-EPI Data

## SENSE Reconstruction

[sense.jl](sense.jl) defines the function `recon`
for reconstructing a 3D volume from SMS-EPI data.


### Getting Started

`recon` depends on two Julia packages:
[MAT.jl](https://github.com/JuliaIO/MAT.jl) and
[EPI3D.jl](https://github.com/HarmonizedMRI/3DEPI/tree/main/recon/sense).
In order to use EPI3D.jl,
you must first download or clone the repository, e.g.,
```bash
$ git clone https://github.com/HarmonizedMRI/3DEPI
```
The packages can be installed
by doing the following from the Julia REPL:
```julia-repl
julia> ] # enter Julia's package prompt
pkg> add MAT

pkg> dev path/to/3DEPI/recon/sense
```
(Make sure to change `path/to/`
to the location where you downloaded or cloned
the 3DEPI repository.)


### Usage

For each new Julia session,
you must first `include` [sense.jl](sense.jl)
to load the necessary packages
and to define the function `recon`:
```julia
include("path/to/SMS-EPI/recon/sense.jl")
```
This needs to be done just once per Julia session.
Then the following will reconstruct an image
given appropriate inputs:
```julia
xhat = recon(data, mask, smap, resolution; outfile)
```
See the documentation for a description of the inputs.


### Documentation

See the docstring in [sense.jl](sense.jl),
or run the following from the Julia REPL
after `include`ing [sense.jl](sense.jl):
```julia-repl
julia> ? # enter Julia's help mode
help?> recon
```
