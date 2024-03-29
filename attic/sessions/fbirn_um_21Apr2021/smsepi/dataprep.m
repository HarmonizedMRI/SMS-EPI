
clear d2d kx2d dcart

% EPI odd/even correction parameters
dly = 0;  % fraction of 4us sample
th0 = 0;     % odd/even dc phase offset

% apply temporal shift (odd/even linear phase correction)
nt = size(kx,1);
kx = interp1(1:nt, kx, (1:nt)+dly);

% reshape 
clear d2d kx2d
for echo = 1:ny   % EPI echo (not dabecho)
	istart = length(gpre) + (echo+3-1)*length(gx1) + length(gpre) + 1;
	istop = istart + length(gx1) - 1;
	d2d(:,:,echo) = dat(istart:istop, :);
	kx2d(:,echo) = kx(istart:istop);
end
d2d = permute(d2d, [1 3 2]);  % [length(gx1) 64 ncoils]

% apply odd/even dc phase offset
for ic = 1:ncoils
	d2d(:,2:2:end,ic) = bsxfun(@times, exp(1i*th0), d2d(:,2:2:end,ic));
end

% interpolate onto cartesian grid along readout
clear dcart;
[~,Ao,dcfo] = reconecho([], nx, [], [], kx2d(:,1), fov); % odd echoes
[~,Ae,dcfe] = reconecho([], nx, [], [], kx2d(:,2), fov); % even echoes
for ic = 1:ncoils
	x = zeros(nx,ny);
	for iy = 1:2:ny
		x(:,iy) = reconecho(d2d(:,iy,ic), nx, Ao, dcfo);
	end
	for iy = 2:2:ny
		x(:,iy) = reconecho(d2d(:,iy,ic), nx, Ae, dcfe);
	end
	dcart(:,:,ic) = fftshift(fft(fftshift(x,1), [], 1),1);
end
dcart(isnan(dcart)) = 0;
dcart = reshape(dcart, [], ncoils);  % [sum(kmask(:)) ncoils]
