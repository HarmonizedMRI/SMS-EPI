% recon mb=6 SMS data with slice-GRAPPA
% 
% In addition to the parameters defined in set_experimental_params.m,
% the following variables must be present in the MATLAB workspace:
%    a           odd/even echo calibration parameters (constant and linear offset)
%    dcal        individual slice k-space, for slice GRAPPA calibration
%    e           struct with exam ids, see getexamids.m
%    And a few others...

nFramesDiscard = 0;

% CAIPI sampling mask
smask = hmriutils.epi.getsamplingmask([1 3 5 1 3 5], nx, etl, mb);
%smask = flipdim(smask, 2);  % for testing negative y (PE) gradient on Siemens

% Initialize slice GRAPPA weights
for p = 1:np
    w{p} = [];
end

% slice order for the np partitions
Z_start = hmriutils.epi.getsliceordering(np);

% loop over frames and reconstruct
Irss = zeros(nx, ny, nz, nFrames);
for ifr = (nFramesDiscard+1):nFrames
    % load raw data for this frame, interpolate to Cartesian grid, 
    % and apply odd/even phase correction
    draw = hmriutils.epi.io.readframe([tmpdir ifn], ifr);
    tic
    dfr = hmriutils.epi.rampsampepi2cart(draw, kxo, kxe, nx, fov(1)*100, 'nufft'); 
    toc
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
        ncalx = 48; ncaly = 48; % setting optimal cal region size is an unsolved problem
        Rx = nx/2-ncalx/2:nx/2+ncalx/2-1;
        Ry = ny/2-ncaly/2:ny/2+ncaly/2-1;
        Ry = Ry - (ny-etl);
        ycal = 0*d_ex;
        ycal(Rx, Ry, :, :) = d_ex(Rx, Ry, :, :);

        % slice GRAPPA recon
        if isempty(w{p})
            K = [5 5];
            [y, w{p}] = hmriutils.epi.slg.recon(ysms, ycal, Z, nz, smask, K);
        else
            y = hmriutils.epi.slg.recon(ysms, ycal, Z, nz, smask, K, 'w', w{p});
        end

        % partial Fourier recon
        %Irss(:,:,Z,ifr) = hmriutils.epi.slg.recon_pfky(y, ny, 'homodyne');
        Irss(:,:,Z,ifr) = hmriutils.epi.slg.recon_pfky(y, ny, 'zerofill');

        % display
        msk = Irss(:,:,:,ifr)>0.0*max(Irss(:));
        if strcmp(E.vendor, 'GE')
            im(flipdim(Irss(:,:,:,ifr).*msk,2)); %, 10*abs(Irss(:,:,:,ifr)-Icalrss).*msk));
        else
            im(Irss(:,:,:,ifr).*msk); %, 10*abs(Irss(:,:,:,ifr)-Icalrss).*msk));
        end
        pth = [E.subject '-' E.site '-' E.scanner '-' E.date '-' E.session];
        title(sprintf('%s, %s, frame %d', pth, ifn, ifr)); pause(0.25);

        % compare with reference image
        %im(cat(1, Icalrss, Irss(:,:,:,ifr)).*msk); %, 10*abs(Irss(:,:,:,ifr)-Icalrss).*msk));
        %title('left: truth; middle: reconstructed; right: 10xdiff'); pause(1);
    end
end

% apply flips to match product sequences
if strcmp(E.vendor, 'Siemens')
    Irss = flipdim(Irss,2);
    %Irss = flipdim(Irss,1);
end
if strcmp(E.vendor, 'GE')
    %Irss = flipdim(Irss,1);
end

% save to .mat file in temp directory
save([tmpdir ifn '.mat'], 'Irss', '-v7.3')

% scale to int16 and write to .nii
m = 1.1*max(Irss(:));
Irss = uint16(round(Irss/m*2^16));
niftiwrite(Irss, [outputdir ifn '.nii']);
