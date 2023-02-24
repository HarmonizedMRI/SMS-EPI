import EPI3D
using MAT: matread, matwrite
include("utils.jl")

"""
    xhat = recon(data, Ω,  smap; outfile, recon_kwargs...)
    xhat = recon(datafile, smap; outfile, recon_kwargs...)

Reconstruct Cartesian SMS data.

# Arguments
- Either
  - `data`: `Array` of size `(ns, nc, np)`
    containing complex multicoil k-space data
  - `Ω`: `BitArray` (mask) of size `(nx, ny, mb)`
    indicating the k-space locations of each of the acquired data samples
- Or
  - `datafile`: `AbstractString` (e.g., `"path/to/data.mat"`)
    indicating a .mat file
    that contains variables `data` and `mask`
    of the forms indicated above
- `smap`: Either
  - `Array` of size `(nx, ny, nz, nc)`
    containing sensitivity maps for each coil,
    or
  - `AbstractString` (e.g., `"path/to/sense.mat"`)
    indicating a .mat file
    that contains the variable `smap`
    of the form indicated above

# Options
- `outfile = "recon.mat"`: `AbstractString` indicating
  where to save the reconstructed image to disk;
  the resulting .mat file will contain one variable, named `xhat`;
  set `outfile = ""` to skip saving the output to disk
- `recon_kwargs...`: Additional keyword arguments passed to `EPI3D.recon`,
  such as:
  - `lambda = 0`: Regularization parameter (no regularization by default)
    (the regularization encourages spatial smoothness)
  - `niter = 50`: Number of iterations of the optimization algorithm to run

# Returns
- `xhat`: `Array` of size `(nx, ny, nz)`
  containing the reconstructed image

## Notes
- `(nx, ny, nz)` is the matrix size of the reconstructed image
- `ns` is the number of acquired k-space samples for each set of
  simultaneously excited slices;
  `ns = count(Ω)`
- `nc` is the number of coils
- `np` is the number of partitions,
  or sets of simultaneously excited slices
- `mb` is the multi-band factor,
  or number of simultaneously excited slices
- `nz == np * mb` must hold
- Each set of simultaneously excited slices is reconstructed
  by treating the data as 3D data,
  hence `ndims(Ω)` is `3`;
  for non-CAIPI sampling schemes,
  `Ω[kx,ky,kz]` is `false` whenever `kz != 1`
- The sampling mask `Ω` is assumed to be the same
  for each set of simultaneously excited slices
- The samples in `data` must be sorted such that
  `d0 = zeros(eltype(data), size(Ω)); d0[Ω] = data[:,1,1]`
  creates the zero-filled k-space matrix
  (for the first coil and first set of simultaneously excited slices),
  and `data[:,1,1] == d0[Ω]`;
  in particular, for non-fly-back EPI acquisitions care is needed
  to reorder the samples in a fly-back manner
  (to be sequential in k-space, not sequential in time)
"""
function recon(data, Ω, smap;
    outfile = "recon.mat",
    lambda = 0,
    recon_kwargs...
)

    # Ensure consistency between the sizes of the inputs
    _ensure_consistency(data, Ω, smap)

    # Allocate memory for the reconstructed image
    xhat = similar(data, size(smap)[1:3]...)

    # Grab the number of partitions, or sets of simultaneously excited slices
    # (e.g., when acquiring 60 slices with a multi-band factor of 2,
    # `npartitions = 30`)
    npartitions = size(data, 3)

    # Reconstruct each set of simultaneously excited slices
    for i = 1:npartitions

        # For the `i`th set of simultaneously excited slices,
        # grab the indexes of the slices as well as the corresponding data and
        # sensitivity maps
        (slices, y, s) = _loop_setup(data, Ω, smap, i)

        # Reconstruct
        (xsms,) = EPI3D.recon(y, s, Ω; λ = lambda, diff_dims = 1:2,
            recon_kwargs...)
        xhat[:,:,slices] = xsms

    end

    # Save the reconstructed image to file
    outfile == "" || matwrite(outfile, Dict("xhat" => xhat))

    return xhat

end

function recon(datafile::AbstractString, args...; kwargs...)

    d = matread(datafile)
    return recon(d["data"], d["mask"], args...; kwargs...)

end

function recon(data, Ω, smapfile::AbstractString; kwargs...)

    return recon(data, Ω, matread(smapfile)["smap"]; kwargs...)

end
