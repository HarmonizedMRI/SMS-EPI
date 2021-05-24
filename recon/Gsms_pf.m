function A = Gsms_pf(KZ, Z, sens, imask, varargin)
% function A = Gsms_pf(KZ, Z, sens, imask, varargin)
%
% SMS EPI system matrix. Assumes 3/4 partial fourier (in ky)
%
% KZ       [ny]            (cycles/cm) kz encoding value along EPI train
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
	'odim', [arg.nx*arg.ny*arg.nc], ...      % size of A*x
	'forw', @A_forw, 'back', @A_back);
	%'mask', true(arg.np,1), ...

return

% x  [arg.np] Non-negative real
% d  [arg.nx*arg.ny*arg.nc] Complex
%
% Data model for each coil:
%  d_ky = P_ky * F * [S(z1)*exp(1i*kz*z1) S(z2)*exp(1i*kz*z2) ...] * exp(1i*th_l) *[x1; x2; ...]
%  d_ky = y(:,ky,c), c = coil, ky = phase encode
%  x1, x2, ..., xmb = slice 1, 2, ..., mb. Mon-negative real.
%  z1, z2, ..., zmb = z location for each slice (cm)
%  kz(ky) = kz encoding (cycles/cm) for the ky^th phase encode
%  th_l = exp(1i*angle(imgl)), where imgl = phase image from central 1/2 of ky-space
%  S: coil sensitivity for coil c
%  F: 2D FFT
%  P_ky: picks out the y^th phase encode
%
function y = A_forw(arg, x)
	x = embed(x, arg.imask);  % [nx ny mb]
	th = 0*x;                 % low-res image phase estimate
	[nx ny nz] = size(x);

	% multiply by low-res phase image
	for iz = 1:size(x,3)
		d = fftshift(fft2(fftshift(x(:,:,iz))));
		d((end/2-nx/2+1):(end/2+nx/2), (end/2-ny/2+1):(end/2+ny/2)) = 0;
		iml = fftshift(ifft2(fftshift(d))); 
		x(:,:,iz) = x .* exp(1i*angle(iml));
	end
	
	% remaining forward operations 
	y = zeros(arg.nx, arg.ny, arg.nc);
	for ic = 1:arg.nc
		for ikzl = 1:length(arg.kzlevels)
			xsum = sum(arg.ekzz(:,:,:,ikzl) .* arg.sens(:,:,:,ic) .* x, 3);
			tmp = fftshift(fft2(fftshift(xsum)));
			y(:,arg.pegroup{ikzl},ic) = tmp(:,arg.pegroup{ikzl}); 
		end
	end
	y = y(:);
return

function x = A_back(arg, y)
	y = reshape(y, [arg.nx arg.ny arg.nc]);
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
	x = real_nonneg(x(arg.imask));  % [arg.np]
return;

function b = real_nonneg(a)
	b = max(real(a), 0);
return

