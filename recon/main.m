% curdir = pwd; cd /opt/matlab/toolbox/irt/; setup; cd(curdir);

%% 4/10/21 commit 69bb5e49ac2b19c19dccc0ffb27e52eb42d187ec (develop branch)

% load data
pfile = 'P30720.7'; 
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [8852 1 20 1 18]
dat = flipdim(dat,1); % yikes

% get readout
addpath ../sequence/
fmri2depi;   % create gpre, gx1, gx. nx = ny = 64; fov = 26; etc
[kx,ky] = toppe.utils.g2k([gx(:) gy(:)]);  % kx = cycles/cm
kx = [kx; zeros(length(gx)-length(kx),1)];

% get data for all echoes in first frame (for which gy is off)
coil = 1;
slice = 3;
frame = 1;

for echo = 1:ny   % EPI echo (not dabecho)
	istart = length(gpre) + (echo-1)*length(gx1) + 1;
	istop = istart + length(gx1) - 1;
	dat1 = dat(istart:istop, coil, slice, 1, frame); 
	kx1 = kx(istart:istop);
	x(:,echo) = reconecho(dat1, kx1, nx, fov, gx1);
end

subplot(121); im(abs(x)); title('mag');
subplot(122); im(angle(x(:,1:2:end)./x(:,2:2:end)), [-pi/3 pi/3]); colorbar;  % paired phase difference

return

% Test 
xtrue = zeros(nx,1);
xtrue(nx/4:(3*nx/4)) = 1;
xtrue(3*nx/8:(5*nx/8)) = 1.5;
d = A * xtrue;               % 'acquired' data
xhat = A' * (d(:) .* dcf(:));   % reconstruct (perform inufft along x)

subplot(131); plot(xtrue); title('true object');
subplot(132); plot(abs(xhat)); title('reconstructed');
subplot(133); plot(angle(xhat)); title('angle');
