% curdir = pwd; cd /opt/matlab/toolbox/irt/; setup; cd(curdir);


% load data
pfile = 'P30720.7';  % 4/10/21 commit 69bb5e49ac2b19c19dccc0ffb27e52eb42d187ec (develop branch)
pfile = 'P31232.7'  % fixed gyblip (factor 2) 4/10/21 commit 64b3fd0ef73bb040957781d400d60487696bb758 (develop branch)
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [8852 1 20 1 18]
dat = flipdim(dat,1); % yikes

% get sequence 
addpath ../sequence/
fmri2depi;   % gpre, gx1, gx. nx = ny = 64; fov = 26; etc
[kx,ky] = toppe.utils.g2k([gx(:) gy(:)]);  % kx = cycles/cm
kx = [kx; zeros(length(gx)-length(kx),1)];

% get data for desired frame/slice and reshape
frame = 10;
slice = 3;
for echo = 1:ny   % EPI echo (not dabecho)
	istart = length(gpre) + (echo-1)*length(gx1) + 1;
	istop = istart + length(gx1) - 1;
	d2d(:,echo) = dat(istart:istop, coil, slice, 1, frame); 
	kx2d(:,echo) = kx(istart:istop);
end

% odd/even phase correction

% recon 2d image
x = reconepi(d2d, kx2d, nx, fov, gx1(:));

return;

% get data for all echoes 
coil = 1;
%frame = 1;  % y gradient off. See ../sequence/fmri2depi.m
%frame = 2;  % y gradient off, x gradient negated

% reconstruct
if false
for frame = 1:2
	for slice = 1:10
		for echo = 1:ny   % EPI echo (not dabecho)
			istart = length(gpre) + (echo-1)*length(gx1) + 1;
			istop = istart + length(gx1) - 1;
			dat1 = dat(istart:istop, coil, slice, 1, frame); 
			kx1 = kx(istart:istop);
			x(:,echo,slice,frame) = reconecho(dat1, kx1, nx, fov, gx1);
		end

		%figure;
		%subplot(121); im(abs(x(:,:,slice))); title(sprintf('mag (slice %d)', slice));
		%subplot(122); im(angle(x(:,1:2:end,slice)./x(:,2:2:end,slice)), [-pi/3 pi/3]); colorbar;  % paired phase difference
	end
end
save x x
else
load x
end
	
% plot
for echo = [1 31]
	for slice = [1 9]
		figure; 
		plot(squeeze(angle(x(:,echo,slice,1)./x(:,echo,slice,2))));
		title(sprintf('echo %d, slice %d', echo, slice));
	end
end

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
