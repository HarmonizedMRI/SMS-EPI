% curdir = pwd; cd /opt/matlab/toolbox/irt/; setup; cd(curdir);

% 4/10/21 commit 69bb5e49ac2b19c19dccc0ffb27e52eb42d187ec (develop branch)
pfile = 'P30720.7'; 
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [8852 1 20 1 18]
dat = flipdim(dat,1); % yikes
addpath ../sequence/
fmri2depi;   % create gpre, gx1, gx. nx = ny = 64; fov = 26; etc

% get data for one echo
coil = 1;
slice = 1;
frame = 1;
echo = 2;   % EPI echo (not dabecho)
istart = length(gpre) + (echo-1)*length(gx1) + 1;
istop = istart + length(gx1) - 1;
dat = dat(istart:istop, coil, slice, 1, frame); 

% kspace
[kx,ky] = toppe.utils.g2k([gx(:) gy(:)]);  % kx = cycles/cm
kx = kx(istart:istop);

% readout gradient trapezoid (with ramp sampling)
% addpath ../sequence/; fmri2depi;
% save gx1 gx1
%load gx1

% calculate kspace
%[kx,~] = toppe.utils.g2k([gx1(:) gx1(:)]);  % kx = cycles/cm
%kx = kx - max(kx)/2;

% System matrix
%nx = 64;
%fov = 22;  % cm
xinit = zeros(nx,1);
nufft_args = {[nx],[6],[2*nx],[nx/2],'minmax:kb'};
mask = true(nx,1);
A = Gmri([fov(1)*kx(:)],mask,'nufft',nufft_args);

% Density compensation function
dcf = gx1(:)/max(gx1);  % simply the gradient (kspace velocity)

% recon (1d)
x = A' * (dat(:) .* dcf(:));   % reconstruct (perform inufft along x)

subplot(121); plot(abs(x)); title('reconstructed');
subplot(122); plot(angle(x)); title('angle');

return;

% Test 
xtrue = zeros(nx,1);
xtrue(nx/4:(3*nx/4)) = 1;
xtrue(3*nx/8:(5*nx/8)) = 1.5;
d = A * xtrue;               % 'acquired' data
xhat = A' * (d(:) .* dcf(:));   % reconstruct (perform inufft along x)

subplot(131); plot(xtrue); title('true object');
subplot(132); plot(abs(xhat)); title('reconstructed');
subplot(133); plot(angle(xhat)); title('angle');
