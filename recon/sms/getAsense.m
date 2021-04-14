function A = getA(arg)
% function A = getA(arg)
%
% Undersampled 3D Cartesian system matrix
%
% arg.imsize   [1 3] 
% arg.imask    [nx ny nz]      image support (logical)
% arg.kmask    [nx ny nz]      logical mask indicating sampled locations
% arg.sens     [nx ny nz nc]  

arg.nc = size(arg.sens,4);

A = fatrix2('arg', arg, ...
	'mask', arg.imask, ... 
	'odim', [arg.imsize arg.nc], ...      % size of A*x
	'forw', @A_forw, 'back', @A_back);

return

% x  [prod(arg.imsize) 1] 
% y  [arg.nc*sum(arg.kmask(:)) 1]
function y = A_forw(arg, x)
	x = reshape(x, arg.imsize);
	y = zeros([arg.imsize arg.nc]);
	for ic = 1:arg.nc
		y(:,:,:,ic) = fftshift(fftn(fftshift(arg.sens(:,:,:,ic).*x)));
	end
	y = reshape(y, arg.nc*sum(arg.kmask(:)), []);
	%y = y(arg.kmask);
return

function x = A_back(arg, y)
	%y = embed(y, arg.kmask);
	x = zeros(arg.imsize);
	for ic = 1:arg.nc
		tmp = conj(arg.sens(:,:,:,ic)).*fftshift(ifftn(fftshift(y(:,:,:,ic))));
		x = x + tmp;
	end
	x = reshape(x, prod(arg.imsize), []); 
return;
