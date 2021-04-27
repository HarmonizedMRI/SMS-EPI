%% Recon using 3d sense nufft, following ../gre3d/main.m

% pfile and frame to reconstruct
pfile = 'P_smsepi.7';
pfile = 'P_epi_mb3.7';
frame = 20; 

% sens maps
load ../gre3d/sens_bart
sens = sens_bart;

% acquisition info:
% npre (acquire samples prior to first echo), nx, ny, fov
if 0
curdir = pwd;
cd tmp
addpath ~/github/pulseq/matlab
smsepi;   % gpre, gx1, gx, gy, nx, ny, fov, etc
cd(curdir);
save scanparams.mat gpre gx1 gx nx ny fov kx seq ex
else
load scanparams.mat
end

% indeces corresponding to EPI train
ntrap = length(gx1);
istart = length(gpre) + 3*ntrap + length(gpre) + 1; % see smsepi.m
istop = istart + ntrap*ny - 1;

% SMS slice info
slSep = ex.sliceSep;
slThick = seq.slThick;
	
% data
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [nfid ncoils nslices 1 ny]
dat = flipdim(dat,1); % as usual
dat = dat(:,:,1,1,frame);  % [nfid ncoils]
[nfid ncoils] = size(dat);

% kspace locations
[~,gx,gy,gz] = toppe.readmod('readout.mod');  % from smsepi_mb3.tar
gamma = 4.2576;  % kHz/G
dt = 4e-3;       % ms
k = dt*gamma*cumsum([gx gy gz]);

% apply gradient delay
ph = [-0.14 -0.8];  % constant (rad) and linear (rad/cycle) odd/even phase offset
k(:,1) = interp1(1:nfid, k(:,1), (1:nfid) + ph(2)/(2*pi));

% crop out dephaser gradients before start of EPI train
dat = dat(istart:istop,:);
k = k(istart:istop,:);

d2d = reshape(dat, [ntrap ny ncoils]);
for iy = 1:2:ny
  % d2d(:,1:2:end,:,isl,:) = exp(+1i*ph(isl,1)/2)*d2d(:,1:2:end,:,isl,:);
  % d2d(:,2:2:end,:,isl,:) = exp(-1i*ph(isl,1)/2)*d2d(:,2:2:end,:,isl,:);
end

% pick out SMS slices from 3D sens map
isl = slSep/slThick;  % should be integer, see smsepi.m
nz = size(sens,3);
IZmb = [nz/2-isl+1, nz/2, nz/2+isl, nz/2+2*isl];
sens = sens(:,:,IZmb,:);
nz = size(sens,3);

nufft_args = {[nx,ny,nz],[6,6,6],[2*nx,2*ny,2*nz],[nx/2,ny/2,nz/2],'minmax:kb'};
mask = true(nx,ny,nz); % Mask for support
%mask(:,:,nz) = false;
L = 6;
%A0 = Gmri([fov*k(:,1) fov*k(:,2) slSep*nz*k(:,3)], ...
A0 = Gmri([fov*k(:,1) fov*k(:,2) fov*k(:,3)], ...
	mask, 'nufft', nufft_args);
A = Asense(A0, sens);

x0 = reshape(A'*dat(:)/(nx*ny*nz), [nx ny nz]);
W = 1; C = 0;
x = qpwls_pcg1(x0, A, W, dat(:), C, ...
                   'niter', 10);
%tic; [x,res] = cgnr_jfn(A, dat(:), x0(:), 15, 1e-7); toc; % Also works
x = reshape(x, [nx ny nz]);
im(x);

return;





%% First attempt, failed
% sensitivity maps
load ../gre3d/sens_bart
sens = sens_bart;
clear sens_bart

ncoils = size(sens,4);

% acquisition info
if 0
curdir = pwd;
cd ../../../sequence/ %tmp
addpath ~/github/pulseq/matlab
smsepi;   % gpre, gx1, gx, gy, nx, ny, fov, etc
cd(curdir);
[kx,ky] = toppe.utils.g2k([gx(:) gy(:)]);  % kx = cycles/cm
kx = [kx; zeros(length(gx)-length(kx),1)];
save scanparams.mat gpre gx1 gx nx ny fov kx seq ex
return;
else
load scanparams.mat
end

% acquired data
pfile = 'P_smsepi.7';
pfile = 'P_epi_mb3.7';
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [nt ncoils nslices 1 nframes]
dat = flipdim(dat,1); % as usual
frame = 20;
dat = dat(:,:,1,1,frame);  % [nfid ncoils]

% matrix size for reconstruction
mb = 4; % multiband/sms factor (number of simultaneous slices)
[nx ny] = size(sens(:,:,1,1));
imsize = [nx ny mb];

% pick out slices from sensitivity map
isl = round(ex.sliceSep/seq.slThick);
nz = size(sens,3);
IZmb = [ nz/2+1-isl, nz/2, nz/2+isl, nz];
sens = sens(:,:,IZmb,:);

% blipped CAIPI sampling pattern
IZ = caipi(ny,3,1);  % NB! if mb=3, kz=4 not sampled (but recon requires mb=even)

% image support
imask = true(imsize);
imask(:,:,end) = false;

% acquired data
dataprep;  % get dcart = [nx*ny ncoils]

% reconstruct
tol = 1e-6;
xhat = reconsms(dcart(:), IZ, imask, sens, tol);
im(xhat);

return
