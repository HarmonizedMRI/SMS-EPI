% function [smap, spgr] = getsense(varargin)
%
% Compute sensitivity maps from SPGR data acquired on a GE scanner.
%
% Options
% -------
% pfile        File location of the P-file containing the SPGR data
%                Default: 'P,b0.7' (because this data can also be used for
%                          B0 field mapping)
% readoutFile  File location of the TOPPE .mod file used for the readout portion
%              of the SPGR scan
%                Default: 'readout.mod'
% slices       Array of slice indices to keep; for example, the SPGR scans used
%              for sensitivity mapping might include a few more slices above
%              and below those acquired during the SMS-EPI scan to avoid slab
%              profile effects and the resulting negative impact on the
%              computed sensitivity maps
%                Default: 1:size(spgr, 3) (keep all slices)
% outfile      String specifying where to save the sensitivity maps;
%              set to '' to skip saving the output to disk
%                Default: 'sense.mat'
%
% Outputs
% -------
% smap         Sensitivity maps
%                Size: [nx, ny, nz, nc]
% spgr         Square root sum-of-squares SPGR image (useful for providing a
%              reference for comparing the reconstructed SMS-EPI images to)
%                Size: [nx, ny, nz]
%
% Note that both smap and spgr are saved to outfile.
function [smap, spgr] = getsense(varargin)

    arg.pfile = 'P,b0.7';
    arg.readoutFile = 'readout.mod';
    arg.slices = [];
    arg.outfile = 'sense.mat';

    arg = toppe.utils.vararg_pair(arg, varargin);

    % Load data from SPGR scan
    [~, spgr, dsense] = toppe.utils.recon3dft(arg.pfile, ...
        'readoutFile', arg.readoutFile, ...
        'echo', 1, ... % Just need the first SPGR image
        'flipFid', false); % Preserve original image orientation
    dsense = dsense(2:4:end,:,:,:); % Discard extra samples in readout

    % Compute sensitivity maps
    fprintf('Computing sensitivity maps (this may take a while)\n');
    tic; sense = bart('ecalib -r 20', dsense); toc;
    smap = bart('slice 4 0', sense); % Equivalently: smap = squeeze(sense(:,:,:,:,1));

    % Remove unneeded slices
    if ~isempty(arg.slices)
        smap = smap(:,:,arg.slices,:);
        spgr = spgr(:,:,arg.slices);
    end

    % Save results
    if ~isempty(arg.outfile)
        save(arg.outfile, 'smap', 'spgr', '-v7.3');
    end

end
