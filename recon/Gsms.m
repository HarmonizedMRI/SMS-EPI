function A = Gsms(KZ, Z, sens, imask)
% function A = Gsms(KZ, Z, sens, imask)
%
% SMS EPI system matrix
%
% KZ       [ny]            (cycles/cm) kz encoding value along EPI train
% Z        [mb]            (cm) slice offsets
%                          mb = number of simultaneous slices (multiband factor)
% sens     [nx ny mb nc]   coil sensitivity maps 
% imask    [nx ny mb]      image support (logical)

arg.KZ = IZ;
arg.Z = Z;
arg.sens = sens;
arg.imask = imask;

[arg.nx arg.ny arg.mb] = size(imask);

arg.nc = size(arg.sens,4);
arg.np = sum(arg.imask(:));    % number of spatial positions (voxels)

arg.nt = arg.nx*arg.ny;  % number of samples per coil

A = fatrix2('arg', arg, ...
	'idim', [arg.np], ...             % size of x
	'odim', [arg.nt*arg.nc], ...      % size of A*x
	'forw', @A_forw, 'back', @A_back);
	%'mask', true(arg.np,1), ...

return

% x  [arg.np] 
% y  [arg.nt*arg.nc]
function y = A_forw(arg, x)
	x = embed(x, arg.imask);  % [nx ny mb]
	y = zeros(arg.nt*arg.nc,1);
	for ic = 1:arg.nc
		xsum = zeros(arg.nx, arg.ny);
		for iz = 1:arg.mb
			xsum = xsum + exp(1i*2*pi*arg.KZ(iz)*arg.Z(iz)) * arg.sens(:,:,iz,ic)) .* x(:,:,iz);
		end
		y2d = fftshift(fftn(fftshift(xsum)));
		y(((ic-1)*arg.nt+1):(ic*arg.nt)) = y2d(:);
	end
return

function x = A_back(arg, y)
	y = reshape(y, [arg.nx arg.ny arg.nc]);
	x = zeros(arg.np,1);
	xfull = fftshift(ifftn(fftshift(tmp3d)))
	for ic = 1:arg.nc
		y1 = y(:,:,ic);
		tmp = zeros(arg.nx, arg.ny, arg.mb);
		for iz = 1:arg.mb
			tmp(:,:,iz) = exp(-1i*2*pi*arg.KZ(iz)*arg.Z(iz)) * conj(arg.sens(:,:,iz,ic)).*fftshift(ifftn(fftshift(y1(:))));
		end
		xfull = xfull + tmp;
	end
	x = xfull(arg.imask);
return;



