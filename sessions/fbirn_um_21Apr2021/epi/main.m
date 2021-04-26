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
cd tmp
addpath ~/github/pulseq/matlab
fmri2depi;   % gpre, gx1, gx, gy, nx, ny, fov, etc
[~,gx,gy] = toppe.readmod('readout.mod');
cd ..
gamma = 4.2576;    % kHz/G
dt = 4e-3;         % ms
kx = gamma*dt*cumsum(gx);
ky = gamma*dt*cumsum(gy);

%save scanparams.mat gpre gx1 gx nx ny fov kx ky
%load scanparams.mat

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

if ~exist('ph', 'var')
ph = getoephase(d2d(:,:,:,:,3:4), kxo, kxe, nx, fov);

figure;
hold on; plot(1:nslices, ph(:,1), 'ro');
plot(1:nslices, ph(:,2), 'go');
plot(1:nslices, ph(:,3), 'bo');
legend('dc', 'x', 'y');
xlabel('slice');
return;
end


%% Apply ph and reconstruct

coil = 20; frame = 8;
% slice = 40;

% First recon the Gmri way
nufft_args = {[ny,nx],[6,6],[2*ny,2*nx],[ny/2,nx/2],'minmax:kb'};
mask = true(ny,nx); % Mask for support
L = 6;
%kx2 = (interp1(1:length(kx), kx, (1:length(kx)) + 0.99))';
A = Gmri([fov*kx((npre+1):(end-1)) fov*ky((npre+1):(end-1))], ...
	mask, 'nufft', nufft_args);
d2ddc = zeros(size(d2d(:,:,coil,slice,frame)));
for iy = 1:2:ny
	d2ddc(:,iy) = d2d(:,iy,coil,slice,frame).*dcfo;
end
for iy = 2:2:ny
	d2ddc(:,iy) = d2d(:,iy,coil,slice,frame).*dcfe;
end

d2ddc = flipdim(d2ddc, 1);
d2ddc = flipdim(d2ddc, 2);
x = reshape(A'*d2ddc(:)/nx, [nx ny]);

% Compare w/ 1d nufft + ift way
x2 = recon2depi(d2d(:,:,coil,slice,frame), kxo, kxe, nx, fov);

% apply ph
ntrap = length(kxo);
k.x = zeros(ntrap, ny);
k.y = zeros(ntrap, ny);
for iy = 1:2:ny
	ktmp = [kxo(1)*ones(2,1); kxo; kxo(end)*ones(2,1)]; % to avoid NaN after interpolation
	tmp = interp1(1:length(ktmp), ktmp, (1:length(ktmp)) - ph(slice,2)/2/(2*pi));
	k.x(:,iy) = tmp(3:(end-2));
	k.y(:,iy) = ones(size(kxo))*ky(npre + ntrap*(iy-1) + round(ntrap/2)) - ph(slice,3)/(2*pi)/2/fov;
	d2ddc(:,iy) = d2ddc(:,iy) * exp(-1i*ph(slice,1)/2);
end
for iy = 2:2:ny
	ktmp = [kxe(1)*ones(2,1); kxe; kxe(end)*ones(2,1)];
	tmp = interp1(1:length(ktmp), ktmp, (1:length(ktmp)) - ph(slice,2)/2/(2*pi));
	k.x(:,iy) = tmp(3:(end-2));
	k.y(:,iy) = ones(size(kxe))*ky(npre + ntrap*(iy-1) + round(ntrap/2)) - ph(slice,3)/(2*pi)/2/fov;
	d2ddc(:,iy) = d2ddc(:,iy) * exp(1i*ph(slice,1)/2);
end
A = Gmri([fov*k.x(:) fov*k.y(:)], ...
	mask, 'nufft', nufft_args);
x3 = reshape(A'*d2ddc(:)/nx, [nx ny]);

im(cat(1, x, x2, x3), [0 0.1*max(abs(x(:)))]); 

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

