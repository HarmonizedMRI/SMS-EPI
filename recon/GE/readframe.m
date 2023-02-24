% function [data, mask] = readframe(frame, mb, varargin)
%
% Read one frame of a (fMRI) SMS-EPI data set acquired on a GE scanner.
%
% Inputs
% ------
% frame        Frame to read
% mb           SMS multiband factor
%
% Options
% -------
% pfile        File location of the P-file containing the SMS-EPI data
%                Default: 'P,fmri_mb6.7' (for mb = 6)
% calFrame     Frame to use for ghost correction
%                Default: 2
% calSlice     Slice to use for ghost correction
%                Default: 5 (really any slice will do)
% epiInfoFile  File location of the .mat file containing the epiInfo struct
%                Default: 'epiInfo.mat'
% outfile      String specifying where to save the read-in frame;
%              set to '' to skip saving the output to disk
%                Default: 'frame012.mat' (for frame = 12)
%
% Outputs
% -------
% data         k-space data for the given frame
%                Size: [nsamples, ncoils, npartitions], where npartitions is
%                      the number of slices in the full reconstructed volume
%                      divided by the multiband factor
% mask         k-space CAIPI sampling mask
%                Size: [nx, ny, mb] Note that this mask is the same for each set
%                      of simultaneously excited slices. Also note that
%                      nsamples = sum(mask(:))
function [data, mask] = readframe(frame, mb, varargin)

    arg.pfile = ['P,fmri_mb' num2str(mb) '.7'];
    arg.calFrame = 2;
    arg.calSlice = 5;
    arg.epiInfoFile = 'epiInfo.mat';
    arg.outfile = sprintf('frame%03d.mat', frame);

    arg = toppe.utils.vararg_pair(arg, varargin);

    %% Ghost correction
    % Get calibration data for ghost correction
    load(arg.epiInfoFile);
    sys = toppe.systemspecs(); % Default settings are okay here
    dcal = hmriutils.epi.loadEPIdata_toppe(arg.pfile, arg.calFrame, epiInfo, sys);

    [nx etl nz nc] = size(dcal); % etl = echo train length (number of y phase encodes)

    xcal = fftshift(ifft(ifftshift(dcal, 1), [], 1), 1); % IFT along readout

    % Get ghost correction parameters for one (arbitrary) slice
    xcal = squeeze(xcal(:,:,arg.calSlice,:)); % [nx, etl, nc]

    % If etl is odd, discard the last echo
    xcal = xcal(:,1:(etl-mod(etl, 2)),:);

    % Get odd/even phase mismatch (linear fit coefficients)
    a = hmriutils.epi.getoephase(xcal, false);

    %% Load and reshape data for later processing
    % Load frame
    d = hmriutils.epi.loadEPIdata_toppe(arg.pfile, frame, epiInfo, sys);

    ny = nx;

    % Apply ghost correction
    dc = hmriutils.epi.epiphasecorrect(d, a);

    % Get CAIPI sampling mask
    mask = false(nx, ny, mb);
    mask(:,epiInfo.mask(:,1:mb)) = true;

    % Zero-pad k-space data to account for CAIPI sampling pattern and partial Fourier
    data = zeros(sum(mask(:)), nc, nz);
    for s = 1:nz
        % Prepare the data for slice s.
        % With no CAIPI sampling scheme, all the k-space data for a slice lies
        % in a single kz plane. However, with CAIPI sampling, each ky line
        % actually comes from a different kz plane. To account for this,
        % first we replicate the k-space data, ending up with mb identical
        % kz planes. Then we use the sampling mask to keep only the k-space
        % samples that were actually acquired in each kz plane.
        % Finally, we flatten the 3D data into a 1D array, keeping only the
        % acquired samples (discarding the zeros). (The mask contains the
        % needed k-spatial information, so data doesn't need to store all the
        % extra zeros.)
        ds = zeros(nx, ny, mb, nc);
        for c = 1:nc
            ds(:,1:etl,:,c) = repmat(dc(:,:,s,c), 1, 1, mb);
        end
        index = 0;
        for kk = 1:mb
            for jj = 1:ny
                for ii = 1:nx
                    if mask(ii,jj,kk)
                        index = index + 1;
                        data(index,:,s) = ds(ii,jj,kk,:);
                    end
                end
            end
        end
    end

    % Save results for later processing
    if arg.outfile ~= ''
        save(arg.outfile, 'data', 'mask', '-v7.3');
    end

end
