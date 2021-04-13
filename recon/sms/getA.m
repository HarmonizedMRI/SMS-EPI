function A = getA(arg)
% SMS EPI system matrix based on 3D k-space encoding (Zahneisen et al MRM 2014).
% Does nufft along readout (ramp sampling), and DFT in ky and kz.
%
% arg.dim   [nx ny nz]   size of image 
% arg.k     [nt nd], nt = number of acquired samples, nd = number of dimensions
% arg.nc    number of receive coils

A = fatrix2('arg', arg, 'odim', [arg.dim arg.nc], ...
	'forw', @A_forw, 'back', @A_back);

return

% x  arg.dim            image
% y  [[arg.dim] arg.nc]   coil data 
function y = A_forw(arg, x)
	% dd
	%x = reshape(x, dim);
	y = fftshift(fftn(fftshift(x)));
	%y = y(:);
return

function x = A_back(arg, y)
	%y = reshape(y, dim);
	x = fftshift(ifftn(fftshift(y)));
	%x = x(:);
return;
