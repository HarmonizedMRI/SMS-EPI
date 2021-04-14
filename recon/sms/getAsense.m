function A = getAsense(arg)
% function A = getAsense(arg)
%
% Undersampled 3D Cartesian system matrix
%
% arg.imsize   [1 3] 
% arg.imask    [nx ny nz]      image support (logical)
% arg.kmask    [nx ny nz]      logical mask indicating sampled locations
% arg.sens     [nx ny nz nc]  

arg.nc = size(arg.sens,4);
arg.nt = sum(arg.kmask(:));    % number of acquired samples (per coil)
arg.np = prod(arg.imsize);     % number of spatial positions (voxels)

A = fatrix2('arg', arg, ...
	'mask', true(arg.np,1), ...
	'odim', [arg.nt*arg.nc], ...      % size of A*x
	'forw', @A_forw, 'back', @A_back);

return

% x  [arg.np 1] 
% y  [arg.nc*arg.nt 1]
function y = A_forw(arg, x)
	x = reshape(x, arg.imsize);
	y = zeros(arg.nt*arg.nc,1);
	for ic = 1:arg.nc
		tmp = fftshift(fftn(fftshift(arg.sens(:,:,:,ic).*x)));
		y(((ic-1)*arg.nt+1):(ic*arg.nt)) = tmp(arg.kmask);
	end
return

function x = A_back(arg, y)
	x = zeros(arg.imsize);
	for ic = 1:arg.nc
		tmp = embed(y(((ic-1)*arg.nt+1):(ic*arg.nt)), arg.kmask);
		tmp = conj(arg.sens(:,:,:,ic)).*fftshift(ifftn(fftshift(tmp)));
		x = x + tmp;
	end
	x = reshape(x, arg.np, []); 
return;
