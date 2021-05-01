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

% precompute exponentials (for speed up)
if isempty(arg.zmap)
	arg.ekzz = zeros(arg.nx, arg.ny, arg.mb, length(arg.kzlevels));
	for ikzl = 1:length(arg.kzlevels)
		for iz = 1:arg.mb
			arg.ekzz(:,:,iz,ikzl) = exp(2i*pi*arg.kzlevels(ikzl)*arg.Z(iz));
		end
	end
else
	arg.ekzzmap = zeros(arg.nx, arg.ny, arg.mb, arg.ny);
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
				xsum = sum(arg.ekzz(:,:,:,ikzl) .* arg.sens(:,:,:,ic) .* x, 3);
				tmp = fftshift(fft2(fftshift(xsum)));
				y(:,arg.pegroup{ikzl},ic) = tmp(:,arg.pegroup{ikzl}); 
			end
		else
			% do one fft for every ky encode (echo)
			for iy = 1:arg.ny
				xsum = sum(arg.ekzzzmap(:,:,:,iy) .* arg.sens(:,:,:,ic) .* x, 3);
				tmp = fftshift(fft2(fftshift(xsum)));
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
				x1 = fftshift(ifft2(fftshift(y1)));
			
				xc = xc + conj(arg.ekzz(:,:,:,ikzl) .* arg.sens(:,:,:,ic)) .* x1;
			end
		else
			for iy = 1:arg.ny
				% P^H
				y1 = zeros(arg.nx, arg.ny);  
				y1(:,iy) = y(:,iy,ic);

				% F^H
				x1 = fftshift(ifft2(fftshift(y1)));
			
				xc = xc + conj(arg.ekzzzmap(:,:,:,iy) .* arg.sens(:,:,:,ic)) .* repmat(x1, [1 1 arg.mb]);
			end
		end
		x = x + xc;
	end
	x = x(arg.imask);  % [arg.np]
return;


% old code

% tried group operations when passing zamp, but it's actually a bit slower

	for ic = 1:arg.nc
		arg.sensrepmat(:,:,:,:,ic) = repmat(arg.sens(:,:,:,ic), [1 1 1 arg.ny]);
	end

	if ~isempty(arg.zmap)
		xrepmat = repmat(x, [1 1 1 arg.ny]);
	end
			for ikzl = 1:length(arg.kzlevels)
				PE = arg.pegroup{ikzl};  % phase-encode indeces for this kz level
				tmp = arg.ekzz(:,:,:,arg.pegroup{ikzl}) .* ...
					repmat(arg.sens(:,:,:,ic), [1 1 1 length(PE)]) .* ...
					repmat(x, [1 1 1 length(PE)]);
		%		tmp = arg.ekzz(:,:,:,PE) .* ...
		%			arg.sensrepmat(:,:,:,PE,ic) .* ...
		%			xrepmat(:,:,:,PE);
				tmp = squeeze(sum(tmp,3));
				tmp = fftshift(fft2(fftshift(tmp)));
				y(:,PE,ic) = tmp(:,PE);
			end


% doing 1d fft followed by one line fft is no faster than 2d fft followed by extracting one row
				%tmp = fftshift(fft(fftshift(xsum), [], 2));    % no faster
				%y(:,iy,ic) = fftshift(fft(fftshift(tmp(:,iy)))); 

