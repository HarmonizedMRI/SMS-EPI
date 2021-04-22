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
