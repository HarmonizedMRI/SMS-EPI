% toy example

% object
dim = [64 64 4];
fov = [24 24 4];   % cm
n = dim(1);
nz = dim(3);
clear xtrue
for iz = 2:(nz-1)
	xtrue(:,:,iz) = phantom(n) * (-1)^(iz+1) * iz/nz;
end
xtrue(n/4:3*n/4,n/4:3*n/4,1) = 1;
xtrue(n/4:3*n/4,n/4:3*n/4,iz+1) = 0.5;

% kspace sampling pattern
deltak = 1./fov;     % cycles/cm
kxrange = (-dim(1)/2:(dim(1)/2-1))*deltak(1);
kyrange = (-dim(2)/2:(dim(2)/2-1))*deltak(2);
kzrange = (-dim(3)/2:(dim(3)/2-1))*deltak(3);
[kx, ky, kz] = ndgrid(kxrange, kyrange, kzrange);

% sensitivity maps
nc = 4;
[x y z] = meshgrid(linspace(-1,1,dim(1)), linspace(-1,1,dim(2)), linspace(-1,1,dim(3))); 
x = x/3;
y = y/3;
z = z/3;
sens(:,:,:,1) = exp(x).*exp(y).*exp(z);
sens(:,:,:,2) = exp(-x).*exp(y).*exp(z);
sens(:,:,:,3) = exp(-x).*exp(y).*exp(-z);
sens(:,:,:,4) = exp(x).*exp(-y).*exp(z);

% system matrix
arg.dim = dim;
arg.k = [kx(:) ky(:) kz(:)];
arg.nc = nc;
arg.sens = sens;
A = getA(arg);

% 'acquired' data
yfull = A*xtrue;
SNR = 2;
yfull = yfull + randn(size(yfull))*mean(abs(yfull(:)))/SNR;

% reconstruct
xhat = A'*yfull;
im(cat(1,xtrue,xhat));
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

