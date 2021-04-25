function ph = getoephase(d2d, kxo, kxe, nx, fov)
% function ph = getoephase(d2d, kxo, kxe, nx, ny, fov)
%
% 2D linear fit to odd/even phase difference for each slice.
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
%  kxo:   [ntrap]   (cycles/cm) sampling locations for odd echoes
%  kxe:   [ntrap]   (cycles/cm) sampling locations for even echoes
%  nx     int       Image size (along x)
%  fov              cm (along x)
% 
% Output:
%  ph     [nslices 3] 
%         ph(:,1)  dc offset (rad)
%         ph(:,2)  x linear phase (cycles/fov)
%         ph(:,3)  y linear phase (cycles/fov)

[ntrap ny ncoils nslices nframes] = size(d2d);

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
		do(:,1:2:end)  = d2d(:,1:2:end,coil,isl,1);
		do(:,2:2:end) = d2d(:,2:2:end,coil,isl,2);
		xo = recon2depi(do, kxo, kxo, nx, fov, Ao, dcfo, Ao, dcfo);

		de = 0*d2d(:,:,1,1,1);
		de(:,1:2:end)  = d2d(:,1:2:end,coil,isl,2);
		de(:,2:2:end) = d2d(:,2:2:end,coil,isl,1);
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

	%figure; im(cat(1,xo,xe));
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

