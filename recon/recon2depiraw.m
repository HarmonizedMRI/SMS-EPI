function recon1(dat, kx, ny, npre, nro, dly, th0)
% Recon ramp-sampled EPI data
% 
% dat   [nt] one EPI readout train
% kx    [nt] cycles/cm
% ny    (int) number of phase encodes
% npre  number of samples in prephaser gradient (before 1st echo)
% nro   number of sample in each readout trapezoid (one EPI echo)
% dly   odd/even delay (fraction of one sample)
% th0   odd/even dc phase offset (radians)

% apply temporal shift (odd/even linear phase correction)
nt = length(kx);
kx = interp1(1:nt, kx, (1:nt)-dly);

% reshape data
d2d = zeros(nro,ny);
j2d = zeros(nro,ny);
for echo = 1:ny   % EPI echo (not dabecho)
	istart = npre + (echo-1)*nro + 1;
	istop = istart + nro - 1;
	d2d(:,echo) = dat(istart:istop);
	kx2d(:,echo) = kx(istart:istop);
end
kxo = kx2d(:,1);
kxe = kx2d(:,2);

% apply odd/even dc phase offset
d2d(:,2:2:end) = bsxfun(@times, exp(1i*th0), d2d(:,2:2:end));

% reconstruct
x = recon2depi(d2d, kxo, kxe, nx, fov);

return;
