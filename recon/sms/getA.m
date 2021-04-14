function A = getA(arg)
% SMS EPI system matrix based on 3D k-space encoding (Zahneisen et al MRM 2014).
% Does nufft along readout (ramp sampling), and DFT in ky and kz.
%
% arg.dim   [nx ny nz]   size of image 
% arg.k     [nt nd], nt = number of acquired samples, nd = number of dimensions
% arg.nc    number of receive coils
% arg.sens  [[arg.dim] arg.nc]  coil sensitivity maps (complex)

A = fatrix2('arg', arg, 'odim', [arg.dim arg.nc], ...
	'imask', true(arg.dim), ...
	'forw', @A_forw, 'back', @A_back);

return

% x  arg.dim              image
% y  [[arg.dim] arg.nc]   coil data 
function y = A_forw(arg, x)
	for ic = 1:arg.nc
		y(:,:,:,ic) = fftshift(fftn(fftshift(arg.sens(:,:,:,ic).*x)));
	end
return

function x = A_back(arg, y)
	x = zeros(arg.dim);
	for ic = 1:arg.nc
		tmp = conj(arg.sens(:,:,:,ic)).*fftshift(ifftn(fftshift(y(:,:,:,ic))));
		x = x + tmp;
	end
return;
