% curdir = pwd; cd /opt/matlab/toolbox/irt/; setup; cd(curdir);

% addpath ../sequence/; fmri2depi;
[kx,~] = toppe.utils.g2k([gx1(:) gx1(:)]);  % kx = cycles/cm
kx = kx - max(kx);

% System matrix
%npix = nx;
%xinit = zeros(npix^2,1);
xinit = zeros(nx,1);
%nufft_args = {[ny,nx],[6,6],[2*ny,2*nx],[ny/2,nx/2],'table',2^12,'minmax:kb'};
nufft_args = {[nx],[6],[2*nx],[nx/2],'minmax:kb'};
mask = true(nx,1); % Mask for support (image space)
L = 6;
A = Gmri([fov(1)*kx(:)],mask,'nufft',nufft_args);
%A = Gmri([fov(1)*kx(:) fov(2)*ky(:)],mask,'nufft',nufft_args);

% Density compensation function
dcf = gx1(2:(end-1))/max(gx1);

% Test 
xtrue = zeros(nx,1);
xtrue(nx/4:(3*nx/4)) = 1;
xtrue(3*nx/8:(5*nx/8)) = 1.5;
d = A * xtrue;               % 'acquired' data
m = A' * (d(:) .* dcf(:));   % reconstruct (perform inufft along x)
plot(abs(m))
