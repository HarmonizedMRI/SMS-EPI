% toy example

% object
imsize = [64 64 6];
fov = [24 24 4];   % cm
n = imsize(1);
nz = imsize(3);
clear xtrue
for iz = 2:(nz-1)
	xtrue(:,:,iz) = phantom(n) * (-1)^(iz+1) * iz/nz;
end
xtrue(n/4:3*n/4,n/4:3*n/4,1) = 1;
xtrue(n/4:3*n/4,n/4:3*n/4,iz+1) = 0.5;

% kspace sampling pattern
deltak = 1./fov;     % cycles/cm
kxrange = (-imsize(1)/2:(imsize(1)/2-1))*deltak(1);
kyrange = (-imsize(2)/2:(imsize(2)/2-1))*deltak(2);
kzrange = (-imsize(3)/2:(imsize(3)/2-1))*deltak(3);
[kx, ky, kz] = ndgrid(kxrange, kyrange, kzrange);

% sensitivity maps
nc = 4;
[x y z] = meshgrid(linspace(-1,1,imsize(1)), linspace(-1,1,imsize(2)), linspace(-1,1,imsize(3))); 
x = x/1;
y = y/1;
z = z/4;
sens(:,:,:,1) = exp(x).*exp(y).*exp(z);
sens(:,:,:,2) = exp(-x).*exp(y).*exp(z);
sens(:,:,:,3) = exp(-x).*exp(y).*exp(-z);
sens(:,:,:,4) = exp(x).*exp(-y).*exp(z);

% 3D cartesian multi-coil system matrix
% Later: add B0 field map
%k = [kx(:) ky(:) kz(:)];
arg.imsize = imsize;
arg.imask = true(imsize);
arg.sens = sens;
arg.kmask = true(imsize);
arg.kmask(:,1:2:end,:) = 0;

arg.nc = size(arg.sens,4);
arg.nt = sum(arg.kmask(:));    % number of acquired samples (per coil)
arg.np = prod(arg.imsize);     % number of spatial positions (voxels)

if 0
	A0 = getA(arg);
	A = Asense(A0, sens);  % need to learn how to use this
else
	A = getAsense(arg);    % explicitly implements forw and back using sens maps
end

% 'acquired' data
yfull = A*xtrue(:);
SNR = 40;
yfull = yfull + randn(size(yfull))*mean(abs(yfull(:)))/SNR;

y = yfull;

if 0   % check coil data
for ic = 1:nc
	tmp = y(((ic-1)*arg.nt+1):(ic*arg.nt));
	tmp = embed(tmp, arg.kmask);
	x = fftshift(ifftn(fftshift(tmp)));
	figure; im(cat(1, xtrue.*sens(:,:,:,ic), x))
end
return;
end

if 0
% check A_back 
xhat = A'*y(:);
xhat = reshape(xhat, [arg.imsize]);
im(cat(1, xtrue(:,:,:,1), xhat/3)); colormap jet;
return
end

% reconstruct
type = 'leak';
nbrs = 4;
chat = 0;
dist_power = 1;
kappa = arg.imask;
%[C, ~] = C2sparse(type, kappa, nbrs, chat, dist_power);

W = diag_sp(ones(size(A,1),1));   % weighting matrix
xinit = zeros(imsize);
%[xhat, info] = qpwls_pcg1(xinit(:), A, W, y, 0, 'niter', 20);
[xhat,res] = cgnr_jfn(A, y, xinit(:), 100);
xhat = reshape(xhat, [arg.imsize]);
im(xhat); colormap jet; 
%figure; im(cat(1, xtrue(:,:,:,1), xhat)); colormap jet; 

return


J = 6; K = 2*dim;
nufft_args = {[dim], [6 6 4], 2*dim, dim/2,'minmax:kb'};
% trick: the system matrix is just the transpose of a SENSE image recon matrix!
A = Gmri_SENSE(k, true(dim), 'fov', fov, 'basis', {'dirac'}, ...
	'nufft', nufft_args, 'sens', sens);

% 'acquired' data
yfull = A*xtrue(:);
SNR = 2;
%yfull = yfull + randn(size(yfull))*mean(abs(yfull(:)))/SNR;

% reconstruct
xhat = A'*yfull(:);
%im(cat(1,xtrue,xhat));
%xinit = zeros(dim);
%[x, info] = qpwls_pcg1(xinit(:), A, diag_sp(ones(size(A,1),1)), yfull(:), lambda*C2, 'niter', 30);

return;

% 'acquired' data
nufft_args = {[n,n,nz],[6,6],[2*n,2*n,2*nz],[n/2,n/2,nz/2],'minmax:kb'};
nufft_args = {dim, [6,6], 2*dim, dim/2, 'minmax:kb'};
mask = true(dim); % Mask for support

nufft_args = {[n,n],[6,6],[2*n,2*n],[n/2,n/2],'minmax:kb'};
mask = true(n,n);
L = 6;
A = Gmri([fov(1)*kx(:) fov(2)*ky(:)],mask,'nufft',nufft_args);
%A = Gmri([fov(1)*kx(:) fov(2)*ky(:) fov(3)*kz(:)], mask, 'nufft', nufft_args);

return;

A = Gmri([fov(1)*kx(:) fov(2)*ky(:) fov(3)*kz(:)], mask, ...
	'L', L, 'nufft', nufft_args);

return;

% Initialize system matrix for each z location
for iz = 1:nz
	if ~isempty(arg.zmap)    % do nufft with fieldmap correction
		if isempty(arg.ti)
			error('You must provide ti (sample times) as well as zmap');
		end

		% Initialize a Gmri object for every slice location
		% This isn't the best implementation because Gmri has a function to
		% update zmap on the fly.
		A{iz} = Gmri([fov(1)*kx(:) fov(2)*ky(:)], mask, ...
            'ti', arg.ti, 'zmap', 1i*2*pi*arg.zmap(:,:,iz), 'L', L, 'nufft', nufft_args);
	end
end


% field-map params
tt = 0:dt:(Nt-1)*dt;tt = tt-(Nt-1)*dt;L = 4;fmap = []; % Hz
% nufft params
J = 6;K = 2*dim;
nufft_args = {[dim dim],[J J],[K K],[dim dim]/2,'minmax:kb'};
gambar = 4257; % gamma/2pi in Hz/g
gam = gambar*2*pi; % gamma in radians/g
% trick: the system matrix is just the transpose of a SENSE image recon matrix!
Gsml = Gmri_SENSE(k,true(dim),'fov',[FOV FOV],'basis',{'dirac'}, ...
	'nufft',nufft_args,'exact',0, ...
	'sens',conj(reshape(sens,[dim*dim Nc]))*(-1i*gam*dt), ...
	'ti',-tt,'L',L,'zmap',2*pi*1i*fmap)';

