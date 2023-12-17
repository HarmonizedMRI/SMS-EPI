% recon mb=6 SMS data with slice-GRAPPA
% 
% In addition to the parameters defined in set_experimental_params.m,
% the following variables must be present in the MATLAB workspace:
%    a           odd/even echo calibration parameters (constant and linear offset)
%    dcal        individual slice k-space, for slice GRAPPA calibration

% file to reconstruct
fn = [datdir pfile_mb6];

% number of frames is determined by the 'runs' parameter on the scanner console
nFrames = 8;

% CAIPI sampling mask
smask = hmriutils.epi.getsamplingmask([1 3 5 1 3 5], nx, ny, mb);

% Initialize slice GRAPPA weights
for p = 1:np
    w{p} = [];
end

% slice order for the np partitions
Z_start = hmriutils.epi.getsliceordering(np);

% loop over frames and reconstruct
Irss = zeros(nx, ny, nz, nFrames);
for ifr = 1:nFrames
    % load raw data for this frame, interpolate to Cartesian grid, 
    % and apply odd/even phase correction
    draw = hmriutils.epi.loadframeraw_ge(fn, etl, np, ifr);   % [nFID etl np nc]
    dfr = hmriutils.epi.rampsampepi2cart(draw, kxo, kxe, nx, fov(1)*100, 'spline'); 
    dfr = hmriutils.epi.epiphasecorrect(dfr, a);    %  [nx etl np nc]

    % loop over partitions (shots, or SMS slice groups)
    for p = 1:length(Z_start) 
        fprintf('Reconstructing partition (slice group) %d of %d\n', p, length(Z_start));

        % slices to recon
        Z = Z_start(p):np:nz;   

        % SMS data for one shot/partition
        ysms = squeeze(dfr(:,:,p,:));   % [nx etl nc]

        % calibration data (acquired without z blips)
        d_ex = dcal(:,:,Z,:);
        ncalx = 48; ncaly = 32; % setting optimal cal region size is an unsolved problem
        Rx = nx/2-ncalx/2:nx/2+ncalx/2-1;
        Ry = ny/2-ncaly/2:ny/2+ncaly/2-1;
        ycal = 0*d_ex;
        ycal(Rx, Ry, :, :) = d_ex(Rx, Ry, :, :);

        % do slice GRAPPA recon
        if isempty(w{p})
            K = [5 5];
            [Irss(:,:,Z,ifr), w{p}] = hmriutils.epi.slg.recon(ysms, ycal, Z, nz, smask, K);
        else
            Irss(:,:,Z,ifr) = hmriutils.epi.slg.recon(ysms, ycal, Z, nz, smask, K, w{p});
        end

        % display
        msk = Icalrss>0.1*max(Icalrss(:));
        im(Irss(:,:,:,ifr).*msk); %, 10*abs(Irss(:,:,:,ifr)-Icalrss).*msk));
        title(sprintf('frame %d', ifr)); pause(0.25);

        % compare with reference image
        %im(cat(1, Icalrss, Irss(:,:,:,ifr)).*msk); %, 10*abs(Irss(:,:,:,ifr)-Icalrss).*msk));
        %title('left: truth; middle: reconstructed; right: 10xdiff'); pause(1);
    end
end

