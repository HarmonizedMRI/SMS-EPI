function A = Gsms_pf(KZ, Z, sens, imask, imlo, varargin)
% function A = Gsms_pf(KZ, Z, sens, imask, imlo, varargin)
%
% SMS EPI system matrix. Assumes 3/4 partial fourier (in ky)
%
% KZ       [ny*pf]         (cycles/cm) kz encoding value along EPI train. pf = Partial Fourier factor (3/4)
% Z        [mb]            (cm) slice offsets
%                          mb = number of simultaneous slices (multiband factor)
% sens     [nx ny mb nc]   coil sensitivity maps 
% imask    [nx ny mb]      image support (logical)
% imlo     [nx ny mb]      imlo = exp(1i*th), where th = low-res image phase
%
% Options:
% zmap     [nx ny mb]      relax_map + 2i*pi*field_map (as in Gmri.m)
%                          included as a exp(-zmap * ti(ky)) term
% ti       [ny]            echo times (at center of each echo in EPI train)
% pf       float           Partial Fourier factor. Default = 3/4
%
% Test function: test_Gsms_pf.m

% default options
arg.zmap = [];
arg.ti = [];
arg.pf = 3/4;

% parse vararg inputs
arg = vararg_pair(arg, varargin);

if ~isempty(arg.zmap) & isempty(arg.ti)
	error(' ''ti'' is required when passing a zmap argument');
end

arg.KZ = KZ;
arg.Z = Z;
arg.sens = sens;
arg.imask = imask;
arg.imlo = imlo;

[arg.nx arg.ny arg.mb] = size(imask);

arg.pfskip = arg.ny*(1-arg.pf);  % number ky encodes to skip (at beginning of echo train)

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
	'odim', [arg.nx*arg.ny*arg.pf*arg.nc], ...      % size of A*x
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
%  th_l = angle(imgl)), where imgl = phase image from central 1/2 of ky-space (assuming PF = 3/4)
%  S: coil sensitivity for coil c
%  F: 2D FFT
%  P_ky: picks out the y^th phase encode
%
function y = A_forw(arg, x)
	x = embed(x, arg.imask);  % [nx ny mb]

	% multiply by low-res phase image
	x = x .* arg.imlo;

	% sensitivity coils and kz encoding
	y = zeros(arg.nx, arg.ny-arg.pfskip, arg.nc);
	for ic = 1:arg.nc
		for ikzl = 1:length(arg.kzlevels)
			xsum = sum(arg.ekzz(:,:,:,ikzl) .* arg.sens(:,:,:,ic) .* x, 3);
			tmp = fftshift(fft2(fftshift(xsum)));
			y(:,arg.pegroup{ikzl},ic) = tmp(:, arg.pfskip + arg.pegroup{ikzl}); 
		end
	end

	% PF sampling
	%if arg.pf < 1
	%	y = y(:, (arg.ny*(1-arg.pf)+1):end, :);
	%end
	
	y = y(:);

return


function x = A_back(arg, y)

	% zero-fill to full size (conjugate of PF sampling)
	tmp = reshape(y, [arg.nx arg.ny-arg.pfskip arg.nc]);
	y = zeros(arg.nx, arg.ny, arg.nc);
	y(:, (arg.pfskip+1):end, :) = tmp;
	%y(:, (end/2-arg.ny/4+1):end, :) = tmp;

	x = zeros(arg.nx, arg.ny, arg.mb);
	for ic = 1:arg.nc
		xc = zeros(arg.nx, arg.ny, arg.mb);
		for ikzl = 1:length(arg.kzlevels)
			% P^H
			y1 = zeros(arg.nx, arg.ny);
			yinds = arg.ny*(1-arg.pf) + arg.pegroup{ikzl};
			y1(:, yinds) = y(:, yinds, ic);

			% F^H
			x1 = fftshift(ifft2(fftshift(y1)));
		
			xc = xc + conj(arg.ekzz(:,:,:,ikzl) .* arg.sens(:,:,:,ic)) .* x1;
		end
		x = x + xc;
	end
	x = x .* conj(arg.imlo);
	x = x(arg.imask);  % [arg.np]

return;

