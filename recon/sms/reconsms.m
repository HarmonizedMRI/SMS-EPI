% toy example
dim = [64 64 4];
n = dim(1);
nz = dim(3);
clear x
for iz = 1:nz
	x(:,:,iz) = phantom(n) * (-1)^(iz+1) * iz/nz;
end
fov = [24 24 4];   % cm

% kspace sampling pattern
deltak = 1./fov;     % cycles/cm
kxrange = (-dim(1)/2:(dim(1)/2-1))*deltak(1);
kyrange = (-dim(2)/2:(dim(2)/2-1))*deltak(2);
kzrange = (-dim(3)/2:(dim(3)/2-1))*deltak(3);
[kx, ky, kz] = ndgrid(kxrange, kyrange, kzrange);

% system matrix
arg.dim = dim(1:2);
arg.mask = true(dim(1:2));
A = getA([kx(:) ky(:)], arg);

% 'acquired' data
x = m(:,:,1);
y = A*x;
xhat = A'*y;
im(cat(1,x,xhat));

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

