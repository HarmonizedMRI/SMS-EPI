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

arg.KZ = KZ;
arg.Z = Z;
arg.sens = sens;
arg.imask = imask;

[arg.nx arg.ny arg.mb] = size(imask);

arg.nc = size(arg.sens,4);
arg.np = sum(arg.imask(:));    % number of spatial positions (voxels)

A = fatrix2('arg', arg, ...
	'idim', [arg.np], ...             % size of x
	'odim', [arg.nx*arg.ny*arg.nc], ...      % size of A*x
	'forw', @A_forw, 'back', @A_back);
	%'mask', true(arg.np,1), ...

return

% x  [arg.np] 
% y  [arg.nx*arg.ny*arg.nc]
function y = A_forw(arg, x)
	x = embed(x, arg.imask);  % [nx ny mb]
	y = zeros(arg.nx, arg.ny, arg.nc);
	for ic = 1:arg.nc
		tmp = zeros(arg.nx, arg.ny);
		for iz = 1:arg.mb
			tmp = tmp + exp(1i*2*pi*arg.KZ(iz)*arg.Z(iz)) * ...
				arg.sens(:,:,iz,ic) .* x(:,:,iz);
		end
		y(:,:,ic) = fftshift(fftn(fftshift(tmp)));
	end
	y = y(:);
return

function x = A_back(arg, y)
	y = reshape(y, [arg.nx arg.ny arg.nc]);
	x = zeros(arg.nx, arg.ny, arg.mb);
	for ic = 1:arg.nc
		tmp = zeros(arg.nx, arg.ny, arg.mb);
		for iz = 1:arg.mb
			tmp(:,:,iz) = exp(-1i*2*pi*arg.KZ(iz)*arg.Z(iz)) * ...
				conj(arg.sens(:,:,iz,ic)) .* fftshift(ifftn(fftshift(y(:,:,ic))));
		end
		x = x + tmp;
	end
	x = x(arg.imask);  % [arg.np]
return;



