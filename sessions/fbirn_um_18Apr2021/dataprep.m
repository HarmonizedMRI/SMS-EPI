
clear d2d kx2d dcart

% acquired data
pfile = 'P_smsepi.7';
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [8292 32 1 1 40] (fid ncoils nslices nechoes nviews)
dat = flipdim(dat,1); % yikes
frame = 10;
dat = dat(:,:,1,1,frame);  % [nfid ncoils]
ncoils = size(dat,2);

% EPI odd/even correction parameters
delay = 0.0; %0*0.16;  % fraction of 4us sample
th0 = 0; %*0.2;     % odd/even dc phase offset

% apply temporal shift (odd/even linear phase correction)
nt = size(dat,1);
for ic = 1:ncoils
	dat(:,ic) = interp1(1:nt, dat(:,ic), (1:nt)+delay);
end

% acquisition info
load tmp/info  % gpre, gx1, fov

% reshape 
for echo = 1:ny   % EPI echo (not dabecho)
	istart = length(gpre) + (echo-1)*length(gx1) + 1;
	istop = istart + length(gx1) - 1;
	d2d(:,:,echo) = dat(istart:istop, :);
	kx2d(:,echo) = kx(istart:istop);
end
kmax = max(kx2d(:,1));
d2d = permute(d2d, [1 3 2]);  % [length(gx1) 64 ncoils]

% apply odd/even dc phase offset
d2d(:,2:2:end) = bsxfun(@times, exp(1i*th0), d2d(:,2:2:end));

% interpolate onto cartesian grid along readout
for ic = 1:ncoils
	for iy = 1:ny
		x(:,iy) = reconecho(d2d(:,iy,ic), kx2d(:,iy), nx, fov, gx1); 
		%x(:,iy) = interp1(kx2d(:,iy), d2d(:,iy,ic), linspace(-kmax,kmax,nx)');
	end
	dcart(:,:,ic) = fftshift(fft(fftshift(x,1), [], 1),1);
	%dcart(:,:,ic) = x;
end
dcart = reshape(dcart, [], ncoils);  % [sum(kmask(:)) ncoils]

