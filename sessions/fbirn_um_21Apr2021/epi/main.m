% curdir = pwd; cd /opt/matlab/toolbox/irt/; setup; cd(curdir);

%% Load 2D EPI data
% load data
pfile = 'P_fmri2depi.7';  
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [nt ncoils nslices 1 nframes]
dat = flipdim(dat,1); % as usual

% Demodulate by observed B0 eddy current, which are mostly linear.
% This should help correct for image-domain spatial shift.
%dat = bsxfun(@times, exp(-1i*b0(:)), dat);

%% get sequence info
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

%% reconstruct
frame = 8;   % part of calibration frames (gx and gy both on, positive)

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

