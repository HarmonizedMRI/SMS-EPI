function x = recon2depi(dat, kxo, kxe, nx, fov)
% 2D EPI reconstruction (singe coil).
% Does nufft along x (ramp sampling), and ift in y.
% Assumes that dat/kx has already undergone odd/even echo correction.
%
% dat:  [nt ny]  data (complex)   
% kxo:  [nt]     (cycles/cm) sampling locations for odd echoes
% kxe:  [nt]     (cycles/cm) sampling locations for even echoes
% nx   image size (int)
% fov  cm (assumes square fov)

[nt ny] = size(dat);

x = zeros(nx,ny);

% inverse nufft along readout
% first get odd/even Gmri matrix
[~,Ao,dcfo] = reconecho([], nx, [], [], kxo(:), fov); % odd echoes
[~,Ae,dcfe] = reconecho([], nx, [], [], kxe(:), fov); % even echoes
for iy = 1:2:ny
	x(:,iy) = reconecho(dat(:,iy), nx, Ao, dcfo);
end
for iy = 2:2:ny
	x(:,iy) = reconecho(dat(:,iy), nx, Ae, dcfe);
end

% ift along pe direction
x = fftshift(fft(fftshift(x,2), [], 2),2);

return

