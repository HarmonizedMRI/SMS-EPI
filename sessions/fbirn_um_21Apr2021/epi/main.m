% curdir = pwd; cd /opt/matlab/toolbox/irt/; setup; cd(curdir);

%% Load 2D EPI data and kspace sampling locations
% d2d = [nfid ny ncoils nslices nframes]
% kxo = [nfid]
% kxe = [nfid]

if ~exist('d2d', 'var')

% load data
pfile = 'P_fmri2depi.7';  
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [nt ncoils nslices 1 nframes]
dat = flipdim(dat,1); % as usual

% reshape into [nt ny ncoils] and get odd/even kspace locations
if 0
cd tmp
addpath ~/github/pulseq/matlab
fmri2depi;   % gpre, gx1, gx, gy, nx, ny, fov, etc
cd ..
[kx,ky] = toppe.utils.g2k([gx(:) gy(:)]);  % kx = cycles/cm
kx = [kx; zeros(length(gx)-length(kx),1)];
save scanparams.mat gpre gx1 gx nx ny fov kx
else
load scanparams.mat
end
npre = length(gpre);
ntrap = length(gx1);
[d2d, kxo, kxe] = gedatreshape(dat, kx, npre, ntrap, ny);

% Demodulate by observed B0 eddy current, which are mostly linear.
% This should help correct for image-domain spatial shift.
%dat = bsxfun(@times, exp(-1i*b0(:)), dat);

% Gmri objects for inverse nufft (ramp sampling)
[~,Ao,dcfo] = reconecho([], nx, [], [], kxo, fov); % odd echoes
[~,Ae,dcfe] = reconecho([], nx, [], [], kxe, fov); % odd echoes

end


%% Estimate odd/even phase difference for each slice
% Calibration data:
% frame(s)       gy    gx
% 1, nframes-4   off   positive
% 2, nframes-3   off   negative
% 3, nframes-2   on    positive
% 4, nframes-1   on    negative
% 5, nframes     off   off

if ~exist('a', 'var')

[X,Y] = ndgrid(((-nx/2+0.5):(nx/2-0.5))/nx, ((-ny/2+0.5):(ny/2-0.5))/ny);

[nfid ny ncoils nslices nframes ] = size(d2d);

% a: 2d plane fit parameters: offset (rad), x linear (cycles/fov), y linear (cycles/fov)
a = zeros(nslices, 3);  
for isl = 1:nslices
	fprintf('Getting odd/even phase difference: slice %d of %d', isl, nslices);
	for ib = 1:60; fprintf('\b'); end;
	th = zeros(nx,ny);
	xsos = zeros(nx,ny);  % sum-of-squares coil combined image (for mask)

	for coil = 1:2:ncoils
		do = 0*d2d(:,:,1,1,1);
		do(:,1:2:end)  = d2d(:,1:2:end,coil,isl,3);
		do(:,2:2:end) = d2d(:,2:2:end,coil,isl,4);
		xo = recon2depi(do, kxo, kxo, nx, fov, Ao, dcfo, Ao, dcfo);

		de = 0*d2d(:,:,1,1,1);
		de(:,1:2:end)  = d2d(:,1:2:end,coil,isl,4);
		de(:,2:2:end) = d2d(:,2:2:end,coil,isl,3);
		xe = recon2depi(de, kxe, kxe, nx, fov, Ae, dcfe, Ae, dcfe);

		xm = (abs(xe) + abs(xo))/2;
		th = th + xm.^2.*exp(1i*angle(xe./xo));

		xsos = xsos + xm.^2;
	end

	th = angle(th);
	xsos = sqrt(xsos);
	mask = xsos > 0.1*max(xsos(:));

	% fit phase difference to 2d plane
	H = [ones(sum(mask(:)),1) X(mask) Y(mask)];  % spatial basis matrix (2d linear)
	a(isl,:) = H\th(mask);  
	%thhat = embed(H*a(isl,:)', mask);
	%figure; im(cat(1, th.*mask, thhat, th.*mask-thhat), 1.0*[-1 1]); colormap hsv; colorbar;
end
fprintf('\n');
save a a

else
load a
end

hold on; plot(1:nslices, a(:,1), 'ro');
plot(1:nslices, a(:,2), 'go');
plot(1:nslices, a(:,3), 'bo');
legend('dc', 'x', 'y');
xlabel('slice');

return;


%% reconstruct
frame = 8;   % part of calibration frames (gx and gy both on, positive)
coil = 10;
slice = 32;
x = recon2depi(d2d(:,:,coil,slice,frame), kxo, kxe, nx, fov);
im(x)

return;

dly = 0.0;
th0 = 0.12;
for coil = 1:size(dat,2)
	fprintf('.');
	for slice = 1:size(dat,3)
		x(:,:,slice,coil) = recon2depiraw(dat(:,coil,slice,1,frame), ...
			kx, [nx ny], fov, length(gpre), length(gx1),	dly, th0); 
	end
end
fprintf('\n');
save x x

return;


%% EPI correction parameters
% plot background signal for different delays and odd/even phase offsets
coil = 10;
slice = 32;
frame = 8;   % part of calibration frames (gx and gy both on, positive)
dly = [0 0]; %-0.2:0.02:0.2;  % fraction of 4us sample
th0 = -0.3:0.02:0.3;  % odd/even dc phase offset
load bg ; % background ROI
clear gnr
for idly = 1:length(dly)
	fprintf('.');
	for ith0 = 1:length(th0)
		x = recon2depiraw(dat(:,coil,slice,1,frame), ...
			kx, [nx ny], fov, length(gpre), length(gx1), ...
			dly(idly), th0(ith0));
		gnr(idly,ith0) = mean(abs(x(bg)));
	end
end
fprintf('\n');
surf(gnr);
return;


%% Get B0 eddy current
if 0
% P32256.7  4/11/2021 commit 926a5a0134c32ca01405a9a058f152b85e85e456
% Body coil
% Excites 7 vertical slices along x, centered around x=0
% See ../sequence/fmri2depi
 
pfile = 'P_fmri2depi.7';  
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [8852 1 20 1 18]
dat = flipdim(dat,1); % as usual

load scanparamsepi.mat  % gpre gx1 gx nx ny fov kx ex

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
end

