function A = Gsms(KZ, Z, sens, imask, varargin)
% function A = Gsms(KZ, Z, sens, imask, varargin)
%
% SMS EPI system matrix
%
% KZ       [ny]            (cycles/cm) kz encoding value along EPI train
%                          Arbitrary, but impacts speed of forw/back operations
% Z        [mb]            (cm) slice offsets
%                          mb = number of simultaneous slices (multiband factor)
% sens     [nx ny mb nc]   coil sensitivity maps 
% imask    [nx ny mb]      image support (logical)
%
% Options:
% zmap     [nx ny mb]      relax_map + 2i*pi*field_map (as in Gmri.m)
%                          included as a exp(-zmap * TE(ky)) term
% ti       [ny]            echo times (msec)
%
% Test function: test_Gsms.m

arg.zmap = [];
arg.ti = [];
arg = vararg_pair(arg, varargin);

if ~isempty(arg.zmap) & isempty(arg.ti)
	error(' ''ti'' is required when passing a zmap argument');
end

arg.KZ = KZ;
arg.Z = Z;
arg.sens = sens;
arg.imask = imask;

[arg.nx arg.ny arg.mb] = size(imask);

arg.nc = size(arg.sens,4);
arg.np = sum(arg.imask(:));    % number of spatial positions (voxels)

% group echoes according to kz-encoding (to speed up forw/back operations)
arg.kzlevels = unique(KZ);
for ikzl = 1:length(arg.kzlevels)
	arg.pegroup{ikzl} = find(KZ == arg.kzlevels(ikzl));
end

A = fatrix2('arg', arg, ...
	'idim', [arg.np], ...             % size of x
	'odim', [arg.nx*arg.ny*arg.nc], ...      % size of A*x
	'forw', @A_forw, 'back', @A_back);
	%'mask', true(arg.np,1), ...

return

% x  [arg.np] 
% d  [arg.nx*arg.ny*arg.nc]
%
% Data model for each coil:
%  d_ky = P_ky * F * [S(z1)*exp(1i*kz*z1) S(z2)*exp(1i*kz*z2) ...] * [x1; x2; ...]
%  d_ky = y(:,ky,c), c = coil, ky = phase encode
%  x1, x2, ..., xmb = slice 1, 2, ..., mb 
%  z1, z2, ..., zmb = z location for each slice (cm)
%  kz(ky) = kz encoding (cycles/cm) for the ky^th phase encode
%  S: coil sensitivity for coil c
%  F: 2D FFT
%  P_ky: picks out the y^th phase encode
%
function y = A_forw(arg, x)
	x = embed(x, arg.imask);  % [nx ny mb]
	y = zeros(arg.nx, arg.ny, arg.nc);
	for ic = 1:arg.nc
		if isempty(arg.zmap)
			% group fft operations by kz encoding level
			for ikzl = 1:length(arg.kzlevels)
				xsum = zeros(arg.nx, arg.ny);
				%xfull = bsxfun(@times, repmat(exp(1i*2*pi*arg.kzlevels(ikzl)*arg.Z), [1 1 arg.mb]), ...
				%	arg.sens(:,:,:,ic) .* x;
				for iz = 1:arg.mb
					xsum = xsum + exp(1i*2*pi*arg.kzlevels(ikzl)*arg.Z(iz)) * ...
					arg.sens(:,:,iz,ic) .* x(:,:,iz);
				end
				tmp = fftshift(fftn(fftshift(xsum)));
				y(:,arg.pegroup{ikzl},ic) = tmp(:,arg.pegroup{ikzl}); 
			end
		else
			% do one fft for every ky encode (echo)
			for iy = 1:arg.ny
				xsum = zeros(arg.nx, arg.ny);
				for iz = 1:arg.mb
					xsum = xsum + exp(1i*2*pi*arg.KZ(iy)*arg.Z(iz)) * ...
					arg.sens(:,:,iz,ic) .* x(:,:,iz) .* exp(-arg.zmap(:,:,iz)*arg.ti(iy));
				end
				tmp = fftshift(fftn(fftshift(xsum)));
				y(:,iy,ic) = tmp(:,iy); 
			end
		end
	end
	y = y(:);
return

function x = A_back(arg, y)
	y = reshape(y, [arg.nx arg.ny arg.nc]);
	x = zeros(arg.nx, arg.ny, arg.mb);
	for ic = 1:arg.nc
		xc = zeros(arg.nx, arg.ny, arg.mb);
		if isempty(arg.zmap)
			for ikzl = 1:length(arg.kzlevels)
				% P^H
				y1 = zeros(arg.nx, arg.ny);
				y1(:,arg.pegroup{ikzl}) = y(:,arg.pegroup{ikzl},ic);

				% F^H
				x1 = fftshift(ifftn(fftshift(y1)));
			
				tmp = zeros(arg.nx, arg.ny, arg.mb);
				for iz = 1:arg.mb
					tmp(:,:,iz) = exp(-1i*2*pi*arg.kzlevels(ikzl)*arg.Z(iz)) * ...
						conj(arg.sens(:,:,iz,ic)) .* x1;
				end
				xc = xc + tmp;
			end
		else
			for iy = 1:arg.ny
				% P^H
				y1 = zeros(arg.nx, arg.ny);
				y1(:,iy) = y(:,iy,ic);

				% F^H
				x1 = fftshift(ifftn(fftshift(y1)));
			
				tmp = zeros(arg.nx, arg.ny, arg.mb);
				for iz = 1:arg.mb
					tmp(:,:,iz) = exp(-1i*2*pi*arg.KZ(iy)*arg.Z(iz)) * ...
						conj(arg.sens(:,:,iz,ic)) .* x1 .* conj(exp(-arg.zmap(:,:,iz)*arg.ti(iy)));
				end
				xc = xc + tmp;
			end
		end
		x = x + xc;
	end
	x = x(arg.imask);  % [arg.np]
return;



