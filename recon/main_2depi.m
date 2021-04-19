% curdir = pwd; cd /opt/matlab/toolbox/irt/; setup; cd(curdir);

%% EPI correction parameters
delay = 0.16;  % fraction of 4us sample
th0 = 0.2;   % odd/even dc phase offset

% frame and slice to reconstruct
rec.frame = 12;
rec.slice = 6;
rec.coil  = 1;

%% Get B0 eddy current
% P32256.7  4/11/2021 commit 926a5a0134c32ca01405a9a058f152b85e85e456
% Body coil
% Excites 7 vertical slices along x, centered around x=0
% See ../sequence/fmri2depi
 
pfile = '../sequence/caldata/P32256.7';  
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [8852 1 20 1 18]
dat = flipdim(dat,1); % yikes

load scanparamsP32256.mat  % gpre gx1 gx nx ny fov kx ex

coil = 1;  % body coil 

frame = 1; % gy off
datn = dat(:,coil,2,1,frame);  % slice at x = -2 cm (ex.spacing = 1cm)
datp = dat(:,coil,6,1,frame);  % slice at x = +2 cm
frame = 5; % gx and gy off
datnref = dat(:,coil,2,1,frame); % slice at x = -2 cm
datpref = dat(:,coil,6,1,frame);
dp = datp./datpref;             % remove B0 field (off resonance) phase
dn = datn./datnref;
b0 = (angle(dp.*dn));


%% Load 2D EPI data and reconstruct
% load data
pfile = 'P30720.7';  % 4/10/21 commit 69bb5e49ac2b19c19dccc0ffb27e52eb42d187ec (develop branch)
pfile = 'P31232.7';  % fixed gyblip (factor 2) 4/10/21 commit 64b3fd0ef73bb040957781d400d60487696bb758 (develop branch)
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [8852 1 20 1 18]
dat = flipdim(dat,1); % yikes

% Demodulate by observed B0 eddy current, which are mostly linear.
% This should help correct for image-domain spatial shift.
%dat = bsxfun(@times, exp(-1i*b0(:)), dat);

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
kx = interp1(1:nt, kx, (1:nt)-delay);

% Get data for desired frame/slice and reshape,
% and odd/even echo calibration data.
for echo = 1:ny   % EPI echo (not dabecho)
	istart = length(gpre) + (echo-1)*length(gx1) + 1;
	istop = istart + length(gx1) - 1;
	d2d(:,echo) = dat(istart:istop, rec.coil, rec.slice, 1, rec.frame);  
	kx2d(:,echo) = kx(istart:istop);
end

% apply odd/even dc phase offset
d2d(:,2:2:end) = bsxfun(@times, exp(1i*th0), d2d(:,2:2:end));

x = reconepi(d2d, kx2d, nx, fov, gx1(:));
im(abs(x));

return;

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
