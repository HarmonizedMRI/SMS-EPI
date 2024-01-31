function imsos = getimsos(d2d, kxo, nx, fov)
% function imsos = getimsos(d2d, kxo, nx, fov)
%
% Get unaliased sum-of-squares image from combining +gx and -gx acquisitions.
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
%  nx     int       Image size (along x)
%  fov              cm (along x)
% 
% Output:
%  imsos  [nx nx nslices]

[ntrap ny ncoils nslices nframes] = size(d2d);

% Gmri objects for inverse nufft (ramp sampling)
[~,Ao,dcfo] = reconecho([], nx, [], [], kxo, fov); % odd echoes

[X,Y] = ndgrid(((-nx/2+0.5):(nx/2-0.5))/nx, ((-ny/2+0.5):(ny/2-0.5))/ny);

mask = false(nx, nx, nslices);

for isl = 1:1:nslices
	fprintf('Getting mask: slice %d of %d', isl, nslices);
	for ib = 1:60; fprintf('\b'); end;

	xsos = zeros(nx,nx); 
	for coil = 1:1:ncoils
		do = 0*d2d(:,:,1,1,1);
		do(:,1:2:end)  = d2d(:,1:2:end,coil,isl,1);
		do(:,2:2:end) = d2d(:,2:2:end,coil,isl,2);
		xo = recon2depi(do, kxo, kxo, nx, fov, Ao, dcfo, Ao, dcfo);
		xsos = xsos + abs(xo).^2;
	end
	imsos(:,:,isl) = sqrt(xsos);
end
fprintf('\n');
