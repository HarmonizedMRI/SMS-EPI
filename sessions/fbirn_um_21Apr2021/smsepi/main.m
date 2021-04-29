% pfile and frame to reconstruct
pfile = 'P_smsepi.7';
pfile = 'P_epi_mb3.7';
frame = 10; 
mb = 3;

% sens maps
load ../gre3d/sens_bart
sens = sens_bart;

% acquisition info:
% gpre, gx1 (readout trap), nx, ny, fov
if 0
curdir = pwd;
cd tmp
addpath ~/github/pulseq/matlab
smsepi;   % gpre, gx1, gx, gy, nx, ny, fov, etc
cd(curdir);
save scanparams.mat gpre gx1 gx nx ny fov kx seq ex
else
load scanparams.mat
end

% SMS slice info
slSep = ex.sliceSep;
slThick = seq.slThick;

% caipi pattern
IZ = caipi(ny,mb,1);
	
% data
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [nfid ncoils nslices 1 ny]
dat = flipdim(dat,1); % as usual
dat = dat(:,:,1,1,frame);  % [nfid ncoils]
[nfid ncoils] = size(dat);

% kspace locations
[~,gx,gy,gz] = toppe.readmod('readout.mod');  % from smsepi_mb3.tar
gamma = 4.2576;  % kHz/G
dt = 4e-3;       % ms
k = dt*gamma*cumsum([gx gy gz]);

% apply gradient delay
ph = [-0.14 -1.0];  % constant (rad) and linear (rad/cycle) odd/even phase offset
k(:,1) = interp1(1:nfid, k(:,1), (1:nfid) + ph(2)/2/(2*pi));

% crop out dephaser gradients before start of EPI train
ntrap = length(gx1);
istart = length(gpre) + 3*ntrap + length(gpre) + 1; % see smsepi.m
istop = istart + ntrap*ny - 1;
dat = dat(istart:istop,:);
k = k(istart:istop,:);

% reshape 
dat = reshape(dat, ntrap, ny, ncoils);

% odd/even echo phase offset
dat(:,1:2:end,:) = exp(+1i*ph(1)/2)*dat(:,1:2:end,:);
dat(:,2:2:end,:) = exp(-1i*ph(1)/2)*dat(:,2:2:end,:);

% pick out SMS slices from 3D sens map
isl = slSep/slThick;  % should be integer, see smsepi.m
nz = size(sens,3);
IZmb = [nz/2-isl, nz/2, nz/2+isl]; %, nz/2+2*isl];
sens = sens(:,:,IZmb+1,:);  % TODO: not sure about the offset here

mb = size(sens,3);

% interpolate onto cartesian grid along readout
clear dcart;
kx2d = reshape(k(:,1), ntrap, ny);
[~,Ao,dcfo] = reconecho([], nx, [], [], kx2d(:,1), fov); % odd echoes
[~,Ae,dcfe] = reconecho([], nx, [], [], kx2d(:,2), fov); % even echoes
for ic = 1:ncoils
	x = zeros(nx,ny);
	for iy = 1:2:ny
		x(:,iy) = reconecho(dat(:,iy,ic), nx, Ao, dcfo);
	end
	for iy = 2:2:ny
		x(:,iy) = reconecho(dat(:,iy,ic), nx, Ae, dcfe);
	end
	dcart(:,:,ic) = fftshift(fft(fftshift(x,1), [], 1),1);
end
dcart(isnan(dcart)) = 0;
dcart = reshape(dcart, [], ncoils);  % [sum(kmask(:)) ncoils]

% kz encoding and slice offsets
kzmax = 0.125; %1/(2*slSep); % cycles/cm
KZ = (IZ-mb/2-0.5)/(mb/2)*kzmax;
KZ = (IZ-mb/2-0.5)/((mb-1)/2)*kzmax;  % TODO: get this from scan design scripts (smsepi.m)
Z = [(-mb/2+0.5):(mb/2-0.5)]*slSep;  % slice locations (cm)

% reconstruct
fprintf('Reconstructing...\n');
imask = true(nx,ny,mb);
A = Gsms(KZ, Z, sens, imask);
xinit = zeros(size(imask));
tol = 1e-6; nitmax = 20;
tic; [xhat,res] = cgnr_jfn(A, dcart(:), xinit(imask), nitmax, tol); toc;
xhat = embed(xhat, imask);
subplot(121); im(xhat)
subplot(122); plot(res, 'o-');

return;







nufft_args = {[nx,ny,mb],[6,6,6],[2*nx,2*ny,2*mb],[nx/2,ny/2,mb/2],'minmax:kb'};
mask = true(nx,ny,mb); % Mask for support
%mask(:,:,mb) = false;
A0 = Gmri(k, mask, ...
	'fov', [fov fov slSep*mb], ...
	'nufft', nufft_args);
A = Asense(A0, sens);

x0 = reshape(A'*dat(:)/(nx*ny*mb), [nx ny mb]);
W = 1; C = 0;
x = qpwls_pcg1(x0, A, W, dat(:), C, ...
                   'niter', 10);
%tic; [x,res] = cgnr_jfn(A, dat(:), x0(:), 15, 1e-7); toc; % Also works
x = reshape(x, [nx ny mb]);
im(x);

return;





%% Synthesize data and recon using 3d sense nufft
if false

% data
pfile = '../epi/P_fmri2depi.7';
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [nfid ncoils nslices 1 nframes]
dat = flipdim(dat,1); % as usual
npredat = 83;

% kspace locations
[~,gx,gy,gz] = toppe.readmod('readout.mod'); 
gamma = 4.2576;  % kHz/G
dt = 4e-3;       % ms
k = dt*gamma*cumsum([gx gy gz]);

nprek = 469;  % > npredat since it does 3 echoes at beginning for epi calibration
ntrap = 101;
nx = 64; ny = 64; fov = 25.6;

[nfid ncoils nslices nframes] = size(squeeze(dat));

% apply gradient delay
ph = [-0.14 -0.8];  % constant (rad) and linear (rad/cycle) odd/even phase offset
nk = size(k,1);
%k(:,1) = interp1((1:nk)', k(:,1), (1:nk)' - ph(2)/(2*pi));

% select frame and MB slices
frame = 8;
IZmb = 13:11:54;  %[15 30 45]; 
d = diff(IZmb);
slThick = 0.4;
slSep = d(1)*slThick;  % cm

dat = dat(:,:,IZmb,1,frame);   % [nfid ncoils mb]
sens = sens(:,:,IZmb,:);

mb = length(IZmb);

% crop to EPI train
istartdat = npredat + 1;
istopdat = istartdat + ntrap*ny - 1;
istartk = nprek + 1;
istopk = istartk + ntrap*ny - 1;
dat = dat(istartdat:istopdat,:,:);
k = k(istartk:istopk,:);

% simulate SMS acquisition
dat = permute(dat, [1 3 2]);  % [ntrap*ny mb ncoils]
d2d = reshape(dat, ntrap, ny, mb, []);  % [ntrap ny mb ncoils]
kzmax = 1/(2*slSep); % cycles/cm
IZ = caipi(ny,mb,1);
kz = zeros(ntrap,ny);
for iy = 1:ny
	kz(:,iy) = kzmax*(IZ(iy)-(mb+1)/2)/(mb/2);
end
for iz = 1:mb
	z = (IZmb(iz) - (nslices/2+1) ) * slThick;
	d2d(:,:,mb,ic) = exp(1i*2*pi*kz*z).*d2d(:,:,mb,ic);
end
dat = squeeze(sum(d2d,3));   % simulated SMS data, [ntrap ny ncoils]
k(:,3) = kz(:);

% apply odd/even phase correction
for iy = 1:2:ny
  % d2d(:,1:2:end,:,isl,:) = exp(+1i*ph(isl,1)/2)*d2d(:,1:2:end,:,isl,:);
  % d2d(:,2:2:end,:,isl,:) = exp(-1i*ph(isl,1)/2)*d2d(:,2:2:end,:,isl,:);
end

% odd echoes only (for testing)
if false
ysamp = 1:2:ny;
dat = dat(:,ysamp,:);
k = reshape(k, [ntrap ny 3]);
k = k(:,ysamp,:);
k = reshape(k, [], 3);
end

% image support
imask = true(nx,ny,mb);

% interpolate onto cartesian grid along readout
clear dcart;
kx2d = reshape(k(:,1), ntrap, ny);
[~,Ao,dcfo] = reconecho([], nx, [], [], kx2d(:,1), fov); % odd echoes
[~,Ae,dcfe] = reconecho([], nx, [], [], kx2d(:,2), fov); % even echoes
for ic = 1:ncoils
	x = zeros(nx,ny);
	for iy = 1:2:ny
		x(:,iy) = reconecho(dat(:,iy,ic), nx, Ao, dcfo);
	end
	for iy = 2:2:ny
		x(:,iy) = reconecho(dat(:,iy,ic), nx, Ae, dcfe);
	end
	dcart(:,:,ic) = fftshift(fft(fftshift(x,1), [], 1),1);
end
dcart(isnan(dcart)) = 0;
dcart = reshape(dcart, [], ncoils);  % [sum(kmask(:)) ncoils]

% reconstruct
tol = 1e-6;
%xhat = reconsms(dcart(:), IZ, imask, sens, tol);
%im(xhat);
%return;

nufft_args = {[nx,ny,mb], [6,6,2], [2*nx,2*ny,2*mb], [nx/2,ny/2,mb/2], 'minmax:kb'};
mask = true(nx,ny,mb); % Mask for support
%ssos = sqrt(sum(abs(sens).^2,4));
%mask(ssos > 1e-3) = false;
A0 = Gmri(k, mask, ...
	'fov', [fov fov slSep*(mb-1)], ...
	'nufft', nufft_args);

A = Asense(A0, sens);

x0 = reshape(A'*dat(:)/(nx*ny*mb), [nx ny mb]);
W = 1; C = 0;
x = qpwls_pcg1(x0(:), A, W, dat(:), C, ...
                   'niter', 10);
%tic; [x,res] = cgnr_jfn(A, dat(:), x0(:), 15, 1e-7); toc; % Also works
x = reshape(x, [nx ny mb]);
im(cat(1,x));

return;

end



%% First attempt with my fatrix based on fft, failed
% sensitivity maps
load ../gre3d/sens_bart
sens = sens_bart;
clear sens_bart

ncoils = size(sens,4);

% acquisition info
if 0
curdir = pwd;
cd ../../../sequence/ %tmp
addpath ~/github/pulseq/matlab
smsepi;   % gpre, gx1, gx, gy, nx, ny, fov, etc
cd(curdir);
[kx,ky] = toppe.utils.g2k([gx(:) gy(:)]);  % kx = cycles/cm
kx = [kx; zeros(length(gx)-length(kx),1)];
save scanparams.mat gpre gx1 gx nx ny fov kx seq ex
return;
else
load scanparams.mat
end

% acquired data
pfile = 'P_smsepi.7';
pfile = 'P_epi_mb3.7';
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [nt ncoils nslices 1 nframes]
dat = flipdim(dat,1); % as usual
frame = 20;
dat = dat(:,:,1,1,frame);  % [nfid ncoils]

% matrix size for reconstruction
mb = 4; % multiband/sms factor (number of simultaneous slices)
[nx ny] = size(sens(:,:,1,1));
imsize = [nx ny mb];

% pick out slices from sensitivity map
isl = round(ex.sliceSep/seq.slThick);
nz = size(sens,3);
IZmb = [ nz/2+1-isl, nz/2, nz/2+isl, nz];
sens = sens(:,:,IZmb,:);

% blipped CAIPI sampling pattern
IZ = caipi(ny,3,1);  % NB! if mb=3, kz=4 not sampled (but recon requires mb=even)

% image support
imask = true(imsize);
imask(:,:,end) = false;

% acquired data
dataprep;  % get dcart = [nx*ny ncoils]

% reconstruct
tol = 1e-6;
xhat = reconsms(dcart(:), IZ, imask, sens, tol);
im(xhat);

return
