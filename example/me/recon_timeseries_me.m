% recon SMS data with slice-GRAPPA
% 
% In addition to the parameters defined in set_experimental_params.m,
% the following variables must be present in the MATLAB workspace:
%    a           odd/even echo calibration parameters (constant and linear offset)
%    dcal        individual slice k-space, for slice GRAPPA calibration
%    e           struct with exam ids, see getexamids.m
%    And a few others...

% nFramesDiscard = 0;

% CAIPI sampling mask
%smask = hmriutils.epi.getsamplingmask([1 3 5 1 3 5], nx, etl, mb);
smask = hmriutils.epi.getsamplingmask(ikz, nx, etl, mb,Ry);
%smask = flipdim(smask, 2);  % for testing negative y (PE) gradient on Siemens

% Initialize slice GRAPPA weights and in-plane GRAPPA weights
for p = 1:np
    w{p} = [];
    w_grappa{p} = [];
end

% slice order for the np partitions
Z_start = hmriutils.epi.getsliceordering(np);

% loop over frames and reconstruct
Irss = zeros(nx, ny, nz, nFrames);
figure;
for ifr = (nFramesDiscard+1):nFrames
    % load raw data for this frame, interpolate to Cartesian grid, 
    % and apply odd/even phase correction
    draw = hmriutils.epi.io.readframe([tmpdir ifn], ifr);
    tic
    dfr = hmriutils.epi.rampsampepi2cart(draw, kxo, kxe, nx, fov(1)*100, 'nufft'); 
    toc
    dfr = hmriutils.epi.epiphasecorrect(dfr, a);    %  [nx etl np nc]
%%
    % loop over partitions (shots, or SMS slice groups)
    for p = 1:length(Z_start) 
        fprintf('Reconstructing partition (slice group) %d of %d\n', p, length(Z_start));

        % slices to recon
        Z = Z_start(p):np:nz;   

        % SMS data for one shot/partition
        ysms = insertZeros(squeeze(dfr(:,(iecho-1)*etl+1:iecho*etl,p,:)),2,Ry-1);   % [nx etl*Ry nc]

        % calibration data (acquired without z blips)
        d_ex = insertZeros(dcal(:,:,Z,:),2,Ry-1);
        ycal = 0*d_ex;
        ycal(Calx, Caly, :, :) = d_ex(Calx, Caly, :, :);

        % slice GRAPPA recon
        if isempty(w{p})
            K = [5 5];
            [y, w{p}] = hmriutils.epi.slg.recon(ysms, ycal, Z, nz, smask, K);
        else
            y = hmriutils.epi.slg.recon(ysms, ycal, Z, nz, smask, K, 'w', w{p});
        end

        % in-plane grappa
        ycal_grappa = dcalgrappa(:,:,Z,:);  
        ycal_grappa = ycal_grappa(Calx_grappa,Caly_grappa,:,:);

        
        y(:,2:Ry:end,:,:) = 0; % non-sampled PE loc might be not stricly zero due to slg; this affects the sampling "mask" in grappa()
        if isempty(w_grappa{p})
            [y,w_grappa{p}] = hmriutils.recon.grappa(y,ycal_grappa,[1,Ry],[3,4]);
        else
            y = hmriutils.recon.grappa(y,ycal_grappa,[1,Ry],[3,4],w_grappa{p});
        end

        % partial Fourier recon
        %Irss(:,:,Z,ifr) = hmriutils.epi.slg.recon_pfky(y, ny, 'homodyne');
        Irss(:,:,Z,ifr) = hmriutils.epi.slg.recon_pfky(y, ny, 'zerofill');

        % display
        msk = Irss(:,:,:,ifr)>0.0*max(Irss(:));
        viewargs = {'row',mb};
        if strcmp(E.vendor, 'GE')
            im(viewargs{:},flipdim(Irss(:,:,:,ifr).*msk,2)); %, 10*abs(Irss(:,:,:,ifr)-Icalrss).*msk));
        else
            im(viewargs{:},Irss(:,:,:,ifr).*msk); %, 10*abs(Irss(:,:,:,ifr)-Icalrss).*msk));
        end
        pth = [E.subject '-' E.site '-' E.scanner '-' E.date '-' E.session];
        title(sprintf('%s, %s, frame %d', pth, ifn, ifr)); pause(0.25);

        % compare with reference image
        %im(cat(1, Icalrss, Irss(:,:,:,ifr)).*msk); %, 10*abs(Irss(:,:,:,ifr)-Icalrss).*msk));
        %title('left: truth; middle: reconstructed; right: 10xdiff'); pause(1);
    end
    %%
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



function B = insertZeros(A, dim, k)
% insertZeros(A, dim, k)
% Inserts k zeros *after each element* along dimension dim.
% Works for any N-D array and keeps trailing zeros.

    if nargin < 3
        k = 1;  % default insert 1 zero
    end

    szA = size(A);
    nd = max(dim, ndims(A));   % ensure dim exists
    szA = [szA, ones(1, nd - numel(szA))];  % pad size with 1s if needed

    n = szA(dim);              % original length along dim
    newLen = n * (k + 1);      % new length after inserting k zeros per element

    % New size (only dimension dim changes)
    szB = szA;
    szB(dim) = newLen;

    % Preallocate B with zeros of same class
    B = zeros(szB, class(A));

    % Fill the positions (every k+1 step)
    idx = repmat({':'}, 1, nd);
    idx{dim} = 1:(k+1):newLen;
    B(idx{:}) = A;
end