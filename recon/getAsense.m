function A = getAsense(kmask, imask, sens)
% function A = getAsense(arg)
%
% Undersampled multi-coil 3D Cartesian system matrix
%
% imask    [nx ny nz]      image support (logical)
% kmask    [nx ny nz]      logical mask indicating sampled locations
% sens     [nx ny nz nc]   coil sensitivity maps 

arg.imsize = size(imask);
arg.kmask = kmask;
arg.imask = imask;
arg.sens = sens;

arg.nc = size(arg.sens,4);
arg.nt = sum(arg.kmask(:));    % number of acquired samples (per coil)
%arg.np = prod(arg.imsize);    % number of spatial positions (voxels)
arg.np = sum(arg.imask(:));    % number of spatial positions (voxels)

A = fatrix2('arg', arg, ...
	'idim', [arg.np], ...      % size of A'*y
	'odim', [arg.nt*arg.nc], ...      % size of A*x
	'forw', @A_forw, 'back', @A_back);
	%'mask', true(arg.np,1), ...

return

% x  [arg.np 1] 
% y  [arg.nc*arg.nt 1]
function y = A_forw(arg, x)
	%x = reshape(x, arg.imsize).*arg.imask;
	x = embed(x, arg.imask);
	y = zeros(arg.nt*arg.nc,1);
	for ic = 1:arg.nc
		tmp = fftshift(fftn(fftshift(arg.sens(:,:,:,ic).*x)));
		y(((ic-1)*arg.nt+1):(ic*arg.nt)) = tmp(arg.kmask);
	end
return

function x = A_back(arg, y)
	x = zeros(arg.np,1);
	for ic = 1:arg.nc
		tmp = embed(y(((ic-1)*arg.nt+1):(ic*arg.nt)), arg.kmask);
		tmp = conj(arg.sens(:,:,:,ic)).*fftshift(ifftn(fftshift(tmp)));
		%x = x + tmp.*arg.imask;
		x = x + tmp(arg.imask);
	end
	%x = reshape(x, arg.np, []); 
return;
