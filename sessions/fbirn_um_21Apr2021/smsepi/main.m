%% sensitivity maps
load ../gre3d/sens_bart
sens = sens_bart;
clear sens_bart

ncoils = size(sens,4);

% matrix size for reconstruction
mb = 4; % multiband/sms factor (number of simultaneous slices)
[nx ny] = size(sens(:,:,1,1));
imsize = [nx ny mb];

% pick out slices from sensitivity map
slSep = 5;  % cm (see smsepi.m)
slThick = 0.2812;
isl = round(slSep/slThick);
nz = size(sens,3);
IZmb = [ nz/2-isl, nz/2, nz/2+isl-1, nz];
sens = sens(:,:,IZmb,:);

% blipped CAIPI sampling pattern
IZ = caipi(ny,3,1);  % NB! Acquisition was mb=3, so kz=4 not sampled (recon requires mb=even)

% image support
imask = true(imsize);
imask(:,:,end) = false;

% acquired data
dataprep;  % dcart = [nx*ny ncoils]

% reconstruct
tol = 1e-6;
xhat = reconsms(dcart(:), IZ, imask, sens, tol);
im(xhat);

return
