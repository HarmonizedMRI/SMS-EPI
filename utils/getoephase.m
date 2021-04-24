function ph = getoephase(pfile, readoutfile, npre, ntrap, nx, ny, fov)
% function ph = getoephase(pfile, readoutfile, npre, ntrap, nx, ny, fov)
%
% 2D linear fit to odd/even phase difference for each slice.
%
% Assumes the following calibration data (see ../sequence/fmri2depi.m):
%  frame(s)       gy    gx
%  1, nframes-4   off   positive
%  2, nframes-3   off   negative
%  3, nframes-2   on    positive
%  4, nframes-1   on    negative
%  5, nframes     off   off
%
% Inputs:
%  pfile         GE pfile
%  readoutfile   .mod file for readout (e.g., 'readout.mod');
%  npre   int     Number of acquired samples before first echo
%  ntrap  int     Number of samples in each single-echo trapezoid
%  nx     int     Image size (along x)
%  ny     int     Number of y phase-encodes 
%  fov    [1]     cm
% 
% Output:
%  pn      [nslices 3] 
%         ph[:,1]  dc offset (rad)
%         ph[:,2]  x linear phase (cycles/fov)
%         ph[:,3]  y linear phase (cycles/fov)

% load data
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [nt ncoils nslices 1 nframes]
dat = flipdim(dat,1); % as usual

% kspace sampling locations
[~,gx] = toppe.readmod(readoutfile);
gamma = 4.2576;  % kHz/G
dt = 4e-3;       % ms
kx = gamma*dt*cumsum(gx);

% reshape into [nt ny ...] and get odd/even kspace locations
[d2d, kxo, kxe] = gedatreshape(dat, kx, npre, ntrap, ny);

[nfid ny ncoils nslices nframes] = size(d2d);

% Gmri objects for inverse nufft (ramp sampling)
[~,Ao,dcfo] = reconecho([], nx, [], [], kxo, fov); % odd echoes
[~,Ae,dcfe] = reconecho([], nx, [], [], kxe, fov); % odd echoes

[X,Y] = ndgrid(((-nx/2+0.5):(nx/2-0.5))/nx, ((-ny/2+0.5):(ny/2-0.5))/ny);

ph = zeros(nslices, 3);  
for isl = 1:1:nslices
	fprintf('Getting odd/even phase difference: slice %d of %d', isl, nslices);
	for ib = 1:60; fprintf('\b'); end;
	th = zeros(nx,ny);
	xsos = zeros(nx,ny);  % sum-of-squares coil combined image (for mask)

	for coil = 1:1:ncoils
		do = 0*d2d(:,:,1,1,1);
		do(:,1:2:end)  = d2d(:,1:2:end,coil,isl,3);
		do(:,2:2:end) = d2d(:,2:2:end,coil,isl,4);
		xo = recon2depi(do, kxo, kxo, nx, fov, Ao, dcfo, Ao, dcfo);

		de = 0*d2d(:,:,1,1,1);
		de(:,1:2:end)  = d2d(:,1:2:end,coil,isl,4);
		de(:,2:2:end) = d2d(:,2:2:end,coil,isl,3);
		xe = recon2depi(de, kxe, kxe, nx, fov, Ae, dcfe, Ae, dcfe);

		xm = (abs(xe) + abs(xo))/2;
		th = th + xm.^2.*exp(1i*angle(xe./xo));

		xsos = xsos + xm.^2;
	end

	th = angle(th);
	xsos = sqrt(xsos);
	mask = xsos > 0.1*max(xsos(:));

	% fit phase difference to 2d plane
	H = [ones(sum(mask(:)),1) X(mask) Y(mask)];  % spatial basis matrix (2d linear)
	ph(isl,:) = H\th(mask);  
end
fprintf('\n');

hold on; plot(1:nslices, ph(:,1), 'ro');
plot(1:nslices, ph(:,2), 'go');
plot(1:nslices, ph(:,3), 'bo');
legend('dc', 'x', 'y');
xlabel('slice');

return;

