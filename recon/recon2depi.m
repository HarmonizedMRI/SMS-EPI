function x = reconepi(dat, kx, nx, fov, gx)
% 2D EPI reconstruction.
% Uses nufft along x (ramp sampling), and ift in y
%
% dat: [nt ny]   data (complex)   
% kx:  [nt ny]   cycles/cm
% nx   image size (int)
% fov  cm (assumes square fov)
% gx   [nt 1]  readout gradient (typically a trapezoid) for one echo in EPI train

[nt ny] = size(dat);

x = zeros(nx,ny);

% inverse nufft along readout
for iy = 1:ny
	x(:,iy) = reconecho(dat(:,iy), kx(:,iy), nx, fov, gx);
end

% ift along pe direction
x = fftshift(fft(fftshift(x,2), [], 2),2);

return

