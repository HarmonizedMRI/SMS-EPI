using MAT: matread
using MIRTjim: jim

"""
    showrecon(reconfile; slices, plot_kwargs...)

Display a plot (heatmap) of reconstructed data.

# Arguments
- `reconfile`: `AbstractString` (e.g., `"path/to/recon.mat"`)
  indicating a .mat file
  that contains a variable `xhat`
  that is the reconstructed image to view

# Options
- `slices = :`: Slices to display (`:` means display all slices);
  this kwarg is ignored for 2D data
- `plot_kwargs...`: Additional keyword arguments passed to `MIRTjim.jim`,
  such as:
  - `title = ""`: Plot title
  - `yflip = true`: Whether to flip the image along y
"""
function showrecon(reconfile; slices = :, plot_kwargs...)

    xhat = matread(reconfile)["xhat"]

    if ndims(xhat) > 2
        xslices = xhat[:,:,slices]
    else
        xslices = xhat
    end

    plt = jim(xslices; plot_kwargs...)
    display(plt)

end
