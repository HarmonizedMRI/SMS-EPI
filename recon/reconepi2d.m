% curdir = pwd; cd /opt/matlab/toolbox/irt/; setup; cd(curdir);

pfile = '
[dat, rdb_hdr] = loadpfile(pfile);

% readout gradient trapezoid (with ramp sampling)
% addpath ../sequence/; fmri2depi;
% save gx1 gx1
load gx1

% calculate kspace
[kx,~] = toppe.utils.g2k([gx1(:) gx1(:)]);  % kx = cycles/cm
kx = kx - max(kx)/2;

% System matrix
nx = 64;
fov = 22;  % cm
xinit = zeros(nx,1);
nufft_args = {[nx],[6],[2*nx],[nx/2],'minmax:kb'};
mask = true(nx,1);
A = Gmri([fov(1)*kx(:)],mask,'nufft',nufft_args);

% Density compensation function
dcf = gx1(1:(end-2))/max(gx1);  % simply the gradient (kspace velocity)

% Test 
xtrue = zeros(nx,1);
xtrue(nx/4:(3*nx/4)) = 1;
xtrue(3*nx/8:(5*nx/8)) = 1.5;
d = A * xtrue;               % 'acquired' data
xhat = A' * (d(:) .* dcf(:));   % reconstruct (perform inufft along x)

subplot(131); plot(xtrue); title('true object');
subplot(132); plot(abs(xhat)); title('reconstructed');
subplot(133); plot(angle(xhat)); title('angle');
