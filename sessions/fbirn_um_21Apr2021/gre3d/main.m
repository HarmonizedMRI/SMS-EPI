pfile = 'P_gre3d.7';

% make sensitivity maps (and later: B0 field map)
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

%% do 3d nufft recon (as a check)

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

nx = ny;   % isotropic resolution and fov
fov = 25.6; % see getparams.m (in gre3d.tar)

nufft_args = {[nx,ny,nz],[6,6,6],[2*nx,2*ny,2*nz],[nx/2,ny/2,nz/2],'minmax:kb'};
mask = true(nx,ny,nz); % Mask for support
L = 6;
A = Gmri([fov*kx(:) fov*ky(:) fov*kz(:)], ...
	mask, 'nufft', nufft_args);
coil = 20;
d1 = dat(:,:,:,coil);
x = reshape(A'*d1(:)/(nx*ny*nz), [nx ny nz]);

im(x);
