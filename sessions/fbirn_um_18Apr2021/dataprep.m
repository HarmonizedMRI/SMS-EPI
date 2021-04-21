
clear d2d kx2d dcart

% acquired data
pfile = 'P_smsepi.7';
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [8292 32 1 1 40] (fid ncoils nslices nechoes nviews)
dat = flipdim(dat,1); % as usual
frame = 10;
dat = dat(:,:,1,1,frame);  % [nfid ncoils]
ncoils = size(dat,2);

% EPI odd/even correction parameters
delay = 0.0; %0*0.16;  % fraction of 4us sample
th0 = 0; %*0.2;     % odd/even dc phase offset

% acquisition info
load tmp/info  % gpre, gx1, fov

% apply temporal shift (odd/even linear phase correction)
nt = size(kx,1);
kx = interp1(1:nt, kx, (1:nt)+delay);

% reshape 
for echo = 1:ny   % EPI echo (not dabecho)
	istart = length(gpre) + (echo-1)*length(gx1) + 1;
	istop = istart + length(gx1) - 1;
	d2d(:,:,echo) = dat(istart:istop, :);
	kx2d(:,echo) = kx(istart:istop);
end
kmax = max(kx2d(:,1));
d2d = permute(d2d, [1 3 2]);  % [length(gx1) 64 ncoils]

% flip echoes
%d2d(:,1:2:end,:) = flipdim(d2d(:,1:2:end,:),1);
%kx2d(:,1:2:end) = flipdim(kx2d(:,1:2:end),1);

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

