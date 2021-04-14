function A = getA(arg)
% function A = getA(arg)
%
% Undersampled 3D Cartesian system matrix
%
% arg.imsize   [1 3] 
% arg.imask    [nx ny nz]      image support (logical)
% arg.kmask    [nx ny nz]      logical mask indicating sampled locations

A = fatrix2('arg', arg, ...
	'mask', arg.imask, ... 
	'odim', arg.imsize, ...      % size of A*x
	'forw', @A_forw, 'back', @A_back);

return

% x  [prod(arg.imsize) 1] 
% y  [sum(arg.kmask(:)) 1]
function y = A_forw(arg, x)
	x = reshape(x, arg.imsize);
	y = fftshift(fftn(fftshift(x)));
	y = y(arg.kmask);
return

function x = A_back(arg, y)
	%y = embed(y, arg.kmask);
	x = fftshift(ifftn(fftshift(y)));
	x = reshape(x, prod(arg.imsize), []); 
return;
