function ph = getoephasecoil(d2d, kxo, kxe, nx, fov, mask)
% function ph = getoephasecoil(d2d, kxo, kxe, nx, fov, mask)
%
% 2D linear fit to odd/even phase difference for each slice and coil.
%
% See also ../sequence/fmri2depi.m which acquires the following calibration data:
%  frame(s)       gy    gx
%  1, nframes-4   off   positive
%  2, nframes-3   off   negative
%  3, nframes-2   on    positive    % 1st frame in d2d
%  4, nframes-1   on    negative    % 2nd frame in d2d
%  5, nframes     off   off
%
% Inputs:
%  d2d    [ntrap ny ncoils nslices 2]
%         frame 1: +gx; frame 2: -gx
%  kxo:   [ntrap]           (cycles/cm) sampling locations for odd echoes
%  kxe:   [ntrap]           (cycles/cm) sampling locations for even echoes
%  nx     int               Image size (along x)
%  fov                      cm (along x)
%  mask   [nx nx nslices]   logical object mask (if not provided, is calculated)
% 
% Output:
%  ph     [nslices 3] 
%         ph(:,1)  dc offset (rad)
%         ph(:,2)  x linear phase (rad/fov)
%         ph(:,3)  y linear phase (rad/fov)

[ntrap ny ncoils nslices nframes] = size(d2d);

% Gmri objects for inverse nufft (ramp sampling)
[~,Ao,dcfo] = reconecho([], nx, [], [], kxo, fov); % odd echoes
[~,Ae,dcfe] = reconecho([], nx, [], [], kxe, fov); % odd echoes

[X,Y] = ndgrid(((-nx/2+0.5):(nx/2-0.5))/nx, ((-ny/2+0.5):(ny/2-0.5))/ny);

ph = zeros(nslices, ncoils, 3);  
for isl = 1:nslices
	fprintf('Getting odd/even phase difference: slice %d of %d', isl, nslices);
	for ib = 1:60; fprintf('\b'); end;

	mask2d = mask(:,:,isl);

	th = zeros(nx,ny);
	xsos = zeros(nx,ny);  % sum-of-squares coil combined image (for mask)

	if(sum(mask2d(:))) > 10
	for coil = 1:1:ncoils
		do = 0*d2d(:,:,1,1,1);
		do(:,1:2:end)  = d2d(:,1:2:end,coil,isl,1);
		do(:,2:2:end) = d2d(:,2:2:end,coil,isl,2);
		xo = recon2depi(do, kxo, kxo, nx, fov, Ao, dcfo, Ao, dcfo);

		de = 0*d2d(:,:,1,1,1);
		de(:,1:2:end)  = d2d(:,1:2:end,coil,isl,2);
		de(:,2:2:end) = d2d(:,2:2:end,coil,isl,1);
		xe = recon2depi(de, kxe, kxe, nx, fov, Ae, dcfe, Ae, dcfe);

		th = angle(xe./xo);

		xm = (abs(xe) + abs(xo))/2;

		% fit phase difference to 2d plane
		H = [ones(sum(mask2d(:)),1) X(mask2d) Y(mask2d)];  % spatial basis matrix (2d linear)
		a = H\th(mask2d);
		ph(isl,coil,:) = a;

		%thhat = embed(H*a(:), mask2d);
		%figure; subplot(211); im(xm);
		%subplot(212); im(cat(1,th, thhat, th-thhat), 1*[-1 1]); colormap hsv;
	end
	end
end
fprintf('\n');

return;

hold on; plot(1:nslices, ph(:,1), 'ro');
plot(1:nslices, ph(:,2), 'go');
plot(1:nslices, ph(:,3), 'bo');
legend('dc', 'x', 'y');
xlabel('slice');
title('odd/even phase difference fit parameters');

return;

