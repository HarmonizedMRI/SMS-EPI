function ph = getoephase1d(d2d, kxo, kxe, nx, fov)
% function ph = getoephase(d2d, kxo, kxe, nx, fov)
%
% Get odd/even phase difference for each slice, from first three echoes
% (without y/z phase-encoding)
%
% See also getoephase.m
%
% Inputs:
%  d2d    [ntrap 3 ncoils nslices]   3 echoes with y/z phase-encodes turned off
%         frame 1: +gx; frame 2: -gx
%  kxo:   [ntrap]   (cycles/cm) sampling locations for odd echoes
%  kxe:   [ntrap]   (cycles/cm) sampling locations for even echoes
%  nx     int       Image size (along x)
%  fov              cm (along x)
% 
% Output:
%  ph     [nslices 3] 
%         ph(:,1)  dc offset (rad)
%         ph(:,2)  x linear phase (rad/fov)

[ntrap ny ncoils nslices] = size(d2d);

% Gmri objects for inverse nufft (ramp sampling)
[~,Ao,dcfo] = reconecho([], nx, [], [], kxo, fov); % odd echoes
[~,Ae,dcfe] = reconecho([], nx, [], [], kxe, fov); % odd echoes

[X,Y] = ndgrid(((-nx/2+0.5):(nx/2-0.5))/nx, ((-ny/2+0.5):(ny/2-0.5))/ny);

ph = zeros(nslices, 2);  
for isl = 1:1:nslices
	fprintf('Getting odd/even phase difference: slice %d of %d', isl, nslices);
	for ib = 1:60; fprintf('\b'); end;

	th = zeros(nx,1);
	xsos = zeros(nx,1);  % sum-of-squares coil combined image (for mask)
	
	for coil = 1:1:ncoils
		x1 = reconecho(d2d(:,1,coil,isl), nx, Ao, dcfo);
		x2 = reconecho(d2d(:,2,coil,isl), nx, Ae, dcfe);
		x3 = reconecho(d2d(:,3,coil,isl), nx, Ao, dcfo);

		tmp = angle(x2./x1.*exp(1i*angle(x1./x3)/2));
		th = th + abs(x1).^2.*exp(1i*tmp);
		xsos = xsos + abs(x1).^2;
	end

	th = angle(th);
	xsos = sqrt(xsos);

	mask = xsos > 0.1*max(xsos(:));

	% fit phase difference to 2d plane
	H = [ones(sum(mask(:)),1) X(mask)];  % spatial basis matrix (2d linear)
	ph(isl,:) = H\th(mask);  

	thhat = embed(H*ph(isl,:)', mask);
	figure; plot([th thhat th-thhat]); legend
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

