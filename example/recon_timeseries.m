% recon mb=6 SMS data with slice-GRAPPA

addpath ~/github/HarmonizedMRI/utils/    % +hmriutils functions

fn = [datdir pfile_mb6];  % smsepi,mb=6,pf=0.80,caipiShiftZ=2.tar
mb = 6;

np = nz/mb;

% Initialize slice GRAPPA weights
for p = 1:np
    w{p} = [];
end

% get odd/even echo calibration parameters
load a

% calibration ('ACS') data, fully sampled 2D EPI
load 2d   % datc

nFrames = 8;
Irss = zeros(nx, ny, nz, nFrames);

% start slice order for the np partitions
Z_start = [1:2:np 2:2:np];
Z_start = [Z_start(1:end-2) Z_start(end) Z_start(end-1)]; % last two shots/partitions are swapped

% loop over frames and reconstruct
for ifr = 8 %:nFrames
    % load raw data
    draw = hmriutils.epi.loadframeraw_ge(fn, etl, np, ifr);   % [nFID etl np nc]
    [nFID etl np nc] = size(draw);   % nFID = number of data samples in ADC window

    % interpolate to Cartesian grid (ramp sampling)
    [kxo, kxe] = getk(sysGE, readoutFile, nFID);   % k-space locations for odd and even echoes
    fprintf('Interpolating... '); tic
    dco = hmriutils.epi.rampsamp2cart(draw(:, 1:2:end, :, :), kxo, nx, fovXcm, 'spline');   % odd echoes
    dce = hmriutils.epi.rampsamp2cart(draw(:, 2:2:end, :, :), kxe, nx, fovXcm, 'spline');   % even echoes
    dfr = zeros(nx, etl, np, nc);
    dfr(:,1:2:end, :, :) = squeeze(dco);
    dfr(:,2:2:end, :, :) = squeeze(dce);
    fprintf(' done (%.2fs)\n', toc);

    % apply odd/even phase correction
    dfr = hmriutils.epi.epiphasecorrect(dfr, a);    %  [nx etl np nc]

    % CAIPI sampling mask
    smask = hmriutils.epi.getsamplingmask([1 3 5 1 3 5], nx, ny, mb);

    % loop over partitions (shots, or SMS slice groups)
    for p = 1:length(Z_start) 
        fprintf('Reconstructing partition (slice group) %d of %d\n', p, length(Z_start));

        % slices to recon
        IZ = Z_start(p):np:nz;   

        % SMS data for one shot/partition
        ysms = squeeze(dfr(:,:,p,:));   % [nx etl nc]

        % calibration data (acquired without z blips)
        d_ex = datc(:,:,IZ,:);
        ncalx = 48; ncaly = 32;
        Rx = nx/2-ncalx/2:nx/2+ncalx/2-1;
        Ry = ny/2-ncaly/2:ny/2+ncaly/2-1;
        ycal = 0*d_ex;
        ycal(Rx, Ry, :, :) = d_ex(Rx, Ry, :, :);

        % get slice GRAPPA weights
        if isempty(w{p})
            K = [5 5];
            [~, w{p}] = hmriutils.epi.slg.recon(ysms, ycal, IZ, nz, smask, K);
        end

        % do slice GRAPPA recon
        Irss(:,:,IZ,ifr) = hmriutils.epi.slg.recon(ysms, ycal, IZ, nz, smask, K, w{p});

        % compare with reference images
        Icalrss = zeros(size(ycal(:,:,:,1)));
        Icalfullrss = zeros(size(ycal(:,:,:,1)));
        for z = 1:mb
            [~, Icalrss(:,:,z)] = toppe.utils.ift3(squeeze(ycal(:,:,z,:)), 'type', '2d');
            [~, Icalfullrss(:,:,z)] = toppe.utils.ift3(squeeze(d_ex(:,:,z,:)), 'type', '2d');
        end
        msk = Icalfullrss>0.1*max(Icalfullrss(:));
        im(cat(1, Icalfullrss, Irss, 10*abs(Irss-Icalfullrss).*msk));
        title('left: truth; middle: reconstructed; right: 10xdiff'); pause(1);
    end
end

return

I = flipdim(flipdim(I,2), 1);
I = abs(I);
I = I/max(I(:));
I = single(I);
save(sprintf('%s/Images.mat', datdir), 'I', '-v7.3');

%load Imb1
%im(cat(1, Imb1, I), [0 0.7])
%tic; espiritreco = bart('pics -l1 -r0.001', du, sens_bart); toc;

