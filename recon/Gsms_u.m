function A = Gsms_u(KZ, IY, Z, sens, imask, varargin)
% function A = Gsms_u(KZ, IY, Z, sens, imask, varargin)
%
% SMS EPI system matrix
% Undersampled in ky: sampling pattern is in IY
%
% KZ       [nysamp]        (cycles/cm) kz encoding values along EPI train
% IY       [nysamp]        ky encoding indices, range is 1:ny
% Z        [mb]            (cm) slice offsets
%                          mb = number of simultaneous slices (multiband factor)
% sens     [nx ny mb nc]   coil sensitivity maps 
% imask    [nx ny mb]      image support (logical)
%
% Options:
% zmap     [nx ny mb]      relax_map + 2i*pi*field_map (as in Gmri.m)
%                          included as a exp(-zmap * ti(ky)) term
% ti       [ny]            echo times (at center of each echo in EPI train)
%
% Test function: test_Gsms_u.m

arg.zmap = [];
arg.ti = [];
arg = vararg_pair(arg, varargin);

if ~isempty(arg.zmap) & isempty(arg.ti)
	error(' ''ti'' is required when passing a zmap argument');
end

arg.KZ = KZ;
arg.IY = IY;
arg.Z = Z;
arg.sens = sens;
arg.imask = imask;

[arg.nx arg.ny arg.mb] = size(imask);

arg.nysamp = length(arg.KZ);
assert(length(arg.IY) == length(arg.KZ));

arg.nc = size(arg.sens,4);
arg.np = sum(arg.imask(:));    % number of spatial positions (voxels)

% group echoes according to kz-encoding (to speed up forw/back operations)
arg.kzlevels = unique(KZ);
for ikzl = 1:length(arg.kzlevels)
	tmp = find(KZ == arg.kzlevels(ikzl));
	arg.pegroup{ikzl} = IY(tmp);
end

% precompute exponentials (for speed up)
if isempty(arg.zmap)
	arg.ekzz = zeros(arg.nx, arg.ny, arg.mb, length(arg.kzlevels));
	for ikzl = 1:length(arg.kzlevels)
		for iz = 1:arg.mb
			arg.ekzz(:,:,iz,ikzl) = exp(2i*pi*arg.kzlevels(ikzl)*arg.Z(iz));
		end
	end
else
	arg.ekzzzmap = zeros(arg.nx, arg.ny, arg.mb, arg.ny);
	for iy = 1:arg.ny
		for iz = 1:arg.mb
			arg.ekzzzmap(:,:,iz,iy) = exp(2i*pi*arg.KZ(iy)*arg.Z(iz))  ...
				.* exp(-arg.zmap(:,:,iz)*arg.ti(iy));
		end
	end
end

A = fatrix2('arg', arg, ...
	'idim', [arg.np], ...             % size of x
	'odim', [arg.nx*arg.nysamp*arg.nc], ...      % size of A*x
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
		% group fft operations by kz encoding level
		for ikzl = 1:length(arg.kzlevels)
			xsum = sum(arg.ekzz(:,:,:,ikzl) .* arg.sens(:,:,:,ic) .* x, 3);
			tmp = fftshift(fft2(fftshift(xsum)));
			y(:, arg.pegroup{ikzl}, ic) = tmp(:, arg.pegroup{ikzl}); 
		end
	end

	% undersample
	y = y(:, arg.IY, :);

	y = y(:);
return

function x = A_back(arg, y)
	% zero-fill to full size
	tmp = reshape(y, [arg.nx arg.nysamp arg.nc]);
	y = zeros(arg.nx, arg.ny, arg.nc);
	y(:, arg.IY, :) = tmp;

	x = zeros(arg.nx, arg.ny, arg.mb);
	for ic = 1:arg.nc
		xc = zeros(arg.nx, arg.ny, arg.mb);
		for ikzl = 1:length(arg.kzlevels)
			% P^H
			y1 = zeros(arg.nx, arg.ny);
			y1(:,arg.pegroup{ikzl}) = y(:,arg.pegroup{ikzl},ic);

			% F^H
			x1 = fftshift(ifft2(fftshift(y1)));
		
			xc = xc + conj(arg.ekzz(:,:,:,ikzl) .* arg.sens(:,:,:,ic)) .* x1;
		end
		x = x + xc;
	end
	x = x(arg.imask);  % [arg.np]
return;

