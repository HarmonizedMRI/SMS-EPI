function A = getAsms(IZ, imask, sens)
%
% SMS EPI system matrix
%
% imask    [nx ny nz]      image support (logical)
% kmask    [nx ny nz]      logical mask indicating sampled locations
% sens     [nx ny nz nc]   coil sensitivity maps 

arg.IZ = IZ;
arg.imask = imask;
arg.sens = sens;

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
% y  [arg.nc*arg.nt]
function y = A_forw(arg, x)
	x = embed(x, arg.imask);
	y = zeros(arg.nt*arg.nc,1);
	y2d = zeros(arg.nx, arg.ny);
	for ic = 1:arg.nc
		tmp = fftshift(fftn(fftshift(arg.sens(:,:,:,ic).*x)));
		for iy = 1:arg.ny
			y2d(:,iy) = tmp(:, iy, arg.IZ(iy));
   	end
		y(((ic-1)*arg.nt+1):(ic*arg.nt)) = y2d(:);
	end
return

function x = A_back(arg, y)
	x = zeros(arg.np,1);
	for ic = 1:arg.nc
		tmp = y(((ic-1)*arg.nt+1):(ic*arg.nt));
		tmp = reshape(tmp, arg.nx, arg.ny);
		tmp3d = zeros(arg.nx, arg.ny, arg.mb);
		for iy = 1:arg.ny
			tmp3d(:,iy,arg.IZ(iy)) = tmp(:, iy);
   	end
		%tmp = embed(y(((ic-1)*arg.nt+1):(ic*arg.nt)), arg.kmask);
		tmp = conj(arg.sens(:,:,:,ic)).*fftshift(ifftn(fftshift(tmp3d)));
		x = x + tmp(arg.imask);
	end
return;
