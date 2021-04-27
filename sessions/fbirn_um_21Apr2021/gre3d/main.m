pfile = 'P_gre3d.7';

%% make sensitivity maps (and later: B0 field map)
if false
system('tar xf gre3d.tar readout.mod');
[ims imsos d] = toppe.utils.recon3dft(pfile, ...   % size(d) = [256 64 64 32] = [nfid ny nz ncoils]
	'readoutFile', 'readout.mod', ...
	'flipfid', true, ...
	'flipim', false, ...
	'echo', 1);

d = d(2:4:end,:,:,:);  % oprbw = 31.25 kHz (decimation = 4)
fprintf('getting bart sens maps...');
tic; sens_bart = bart('ecalib -r 20', d); toc;   % takes 14 min
fprintf('\n');
sens_bart = sens_bart(:,:,:,:,1);
save sens_bart sens_bart
end

%% do 3d sense nufft recon (for testing/learning)

% sens maps
% flip to match nufft (Gmri) orientation
load sens_bart;
sens = sens_bart;
sens = flipdim(sens,2);
sens = flipdim(sens,3);

% data
[~,gx,gy,gz,desc,hdrint] = toppe.readmod('readout.mod');  % from gre3d.tar
npre = hdrint(1);
nro = hdrint(2);   % number of 4us samples on gradient plateau
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [nfid ncoils nslices 1 ny]
dat = flipdim(dat,1); % as usual
dat = dat((npre+1):(npre+nro), :, :, :, :);
dat = permute(dat, [1 5 3 2 4]);  % [nro ny nz ncoils]

[nro ny nz ncoils] = size(dat);

% kspace locations
gamma = 4.2576;  % kHz/G
dt = 4e-3;       % ms
k = dt*gamma*cumsum([gx gy gz]);
k = k((npre+1):(npre+nro),:);

kx = repmat(k(:,1), [1 ny nz]);

ky = repmat(k(:,2), [1 ny]);
ky = bsxfun(@times, 2/ny*[(-ny/2+0.5):(ny/2-0.5)], ky);
ky = repmat(ky, [1 1 nz]);

kz = repmat(k(:,3), [1 nz]);
kz = bsxfun(@times, 2/nz*[(-nz/2+0.5):(nz/2-0.5)], kz);
kz = repmat(kz, [1 1 ny]);
kz = permute(kz, [1 3 2]);

% reduce size (ok since delta_k is 1/4 of design fov)
dat = dat(2:4:end,:,:,:);
kx = kx(2:4:end,:,:,:);
ky = ky(2:4:end,:,:,:);
kz = kz(2:4:end,:,:,:);

% reduce matrix size to speed up this test
n = nz/2;
r = (n/2+1):(3*n/2);
sensorig = sens;
sens = zeros(n,n,n,ncoils);
for ic = 1:ncoils
	tmp = fftshift(ifftn(fftshift(sensorig(:,:,:,ic))));
	sens(:,:,:,ic) = fftshift(fftn(fftshift(tmp(r,r,r))));
end
dat = dat(r,r,r,:);
kx = kx(r,r,r);
ky = ky(r,r,r);
kz = kz(r,r,r);

% undersample
ysamp = 1:3:n;
zsamp = 1:3:n;
dat = dat(:,ysamp,zsamp,:);
kx = kx(:,ysamp,zsamp);
ky = ky(:,ysamp,zsamp);
kz = kz(:,ysamp,zsamp);

fov = 25.6; % see getparams.m (in gre3d.tar)

nufft_args = {[n,n,n],[6,6,6],[2*n,2*n,2*n],[n/2,n/2,n/2],'minmax:kb'};
mask = true(n,n,n); % Mask for support
L = 6;
A0 = Gmri([fov*kx(:) fov*ky(:) fov*kz(:)], ...
	mask, 'nufft', nufft_args);
A = Asense(A0, sens);

x0 = reshape(A'*dat(:)/(n*n*n), [n n n]);
W = 1; C = 0;
x = qpwls_pcg1(x0, A, W, dat(:), C, ...
                   'niter', 10);
%tic; [x,res] = cgnr_jfn(A, dat(:), x0(:), 15, 1e-7); toc; % Also works
x = reshape(x, [n n n]);
im(x);

% do single-coil recon as a check
%coil = 20;
%d1 = dat(:,:,:,coil);
%x1 = reshape(A0'*d1(:)/(n*n*n), [n n n]);
%im(x1)
