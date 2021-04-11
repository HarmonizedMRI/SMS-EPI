% curdir = pwd; cd /opt/matlab/toolbox/irt/; setup; cd(curdir);

%% Load eddy current calibration scan
% P32256.7  4/11/2021 commit 
% Excite 7 vertical slices along x. Analyze frames [1 2 nframes-3 nframes] (gy off)
% See ../sequence/fmri2depi
% 
pfile = '../sequence/caldata/P32256.7';  
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [8852 1 20 1 18]
dat = flipdim(dat,1); % yikes

load scanparamsP31744.mat  % gpre gx1 gx nx ny fov kx. Slice separaton = 1 cm.

coil = 1;
frame = 1;
slice = 2;
datp = dat(:,coil,3,1,frame);
datn = dat(:,coil,5,1,frame);
b0 = unwrap(angle(datp.*datn));

return;

%% Load 2D EPI data and reconstruct
% load data
pfile = 'P30720.7';  % 4/10/21 commit 69bb5e49ac2b19c19dccc0ffb27e52eb42d187ec (develop branch)
pfile = 'P31232.7';  % fixed gyblip (factor 2) 4/10/21 commit 64b3fd0ef73bb040957781d400d60487696bb758 (develop branch)
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [8852 1 20 1 18]
dat = flipdim(dat,1); % yikes

% get sequence 
if 0
addpath ../sequence/
fmri2depi;   % gpre, gx1, gx. nx = ny = 64; fov = 26; etc
[kx,ky] = toppe.utils.g2k([gx(:) gy(:)]);  % kx = cycles/cm
kx = [kx; zeros(length(gx)-length(kx),1)];
save scanparamsP31232.mat gpre gx1 gx nx ny fov kx
else
load scanparamsP31232.mat
end

% apply temporal shift (odd/even linear phase correction)
nt = length(kx);
kx = interp1(1:nt, kx, (1:nt)-0.2);

% Get data for desired frame/slice and reshape,
% and odd/even echo calibration data.
frame = 10;
slice = 3;
coil = 1;
for echo = 1:ny   % EPI echo (not dabecho)
	istart = length(gpre) + (echo-1)*length(gx1) + 1;
	istop = istart + length(gx1) - 1;
	d2d(:,echo) = dat(istart:istop, coil, slice, 1, frame);  
	kx2d(:,echo) = kx(istart:istop);

	drefpos = dat(istart:istop, coil, slice, 1, 1);  % y gradient off 
	drefneg = dat(istart:istop, coil, slice, 1, 2);  % y gradient off, gx negated
	xrefpos(:,echo) = reconecho(drefpos, kx2d(:,echo), nx, fov, gx1(:));  % 1d echo profile
	xrefneg(:,echo) = reconecho(drefneg, kx2d(:,echo), nx, fov, gx1(:));
end

% odd/even constant phase offset correction 
% linear fit to pos/neg phase difference, then apply corresponding kspace shift
%te = 4e-6 * length(gx1);  % sec
%f = angle(x(:,end)./x(:,1))/(2*pi*te*(ny-1));  % B0 field map, Hz
%th = angle(xrefpos./xrefneg)/2;
for echo = 1:1:ny   % EPI echo (not dabecho)
	th = angle(xrefpos(:,echo)./xrefneg(:,echo));
	mask = abs(xrefpos(:,echo)) > 0.25*max(abs(xrefpos(:,echo)));
	b0(echo) = mean(th(mask));
	d2d(:,echo) = exp((-1)^(echo)*1i*b0(echo))*d2d(:,echo);
end
for echo = 2:2:ny   % EPI echo (not dabecho)
	%th = angle(xrefpos(:,echo)./xrefneg(:,echo));
	%mask = abs(xrefpos(:,echo)) > 0.2*max(abs(xrefpos(:,echo)));
	%b0(echo) = mean(th(mask));
	%d2d(:,echo) = exp(1i*b0(echo))*d2d(:,echo);
end
%b0 = mean(b0(1:2:end));
%d2d = exp(-1i*b0)*d2d;

% recon 2d image
x = reconepi(d2d, kx2d, nx, fov, gx1(:));
im(abs(x));

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
