# size(data) = (nsamples, ncoils, npartitions)
# size(Ω) = (nx, ny, mb)
# size(smap) = (nx, ny, nz, ncoils)
function _ensure_consistency(data, Ω, smap)

    size(data, 1) == count(Ω) || error("the k-space data contains \
        `size(data, 1)` = $(size(data, 1)) samples, but the sampling mask \
        only specifies sampling locations for `count(Ω)` = $(count(Ω)) samples")

    size(data, 2) == size(smap, 4) || error("the number of coils in the data \
        (`size(data, 2)` = $(size(data, 2)))
        does not match the number of coils in the sensitivity maps \
        (`size(smap, 4)` = $(size(smap, 4))")

    size(data, 3) * size(Ω, 3) == size(smap, 3) || error("the number of \
        reconstructed slices (`size(smap, 3)` = $(size(smap, 3))) \
        must equal the multi-band factor (`size(Ω, 3)` = $(size(Ω, 3))) \
        times the number of sets of simultaneously excited slices \
        (`size(data, 3)` = $(size(data, 3)))")

    all(size(Ω)[1:2] .== size(smap)[1:2]) || error("nx and ny not consistent \
        between the sampling mask $(size(Ω)[1:2]) and the sensitivity maps \
        $(size(smap)[1:2])")

end

# `i` is the loop index, ranging from `1` to `size(data, 3)`
function _loop_setup(data, Ω, smap, i)

    # Grab the total number of reconstructed slices and coils
    (_, _, nz, nc) = size(smap)

    # Specify what slices to reconstruct
    # This assumes the simultaneously excited slices are equally spaced
    slices = i:size(data, 3):nz

    # The expected multi-band factor should be correct
    @assert length(slices) == size(Ω, 3)

    # Grab the data for the set of slices to reconstruct and reshape into a
    # format expected by `EPI3D.recon`
    y = [data[:,c,i] for c = 1:nc]

    # Perform phase correction to account for the fact that the simultaneously
    # excited slices were not necessarily at isocenter
    EPI3D.phasecorrect!(y, i, nz, Ω)

    # Grab the sensitivity maps for the set of slices to reconstruct and reshape
    # into a format expected by `EPI3D.recon`
    s = [smap[:,:,slices,c] for c = 1:nc]

    return (slices, y, s)

end
