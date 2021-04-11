function x = reconepi(dat, kx, nx, fov, gx)
% 2D EPI reconstruction.
% Uses nufft along x (ramp sampling), and ift in y
%
% dat: [nt ny]   
% kx:  [nt 1]  (cycles/cm)
% nx   image size (int)
% fov  cm (assumes square fov)
% gx   [nt 1]  readout gradient (typically a trapezoid) for one echo in EPI train

[nt ny] = size(dat);

x = zeros(nx,ny);

% inufft along readout
for iy = 1:ny
	x(:,iy) = reconecho(dat(:,iy), kx(:,iy), nx, fov, gx);
end

% ift along pe direction
x = fftshift(fft(fftshift(x,2), [], 2),2);

return

% Test 
xtrue = zeros(nx,1);
xtrue(nx/4:(3*nx/4)) = 1;
xtrue(3*nx/8:(5*nx/8)) = 1.5;
d = A * xtrue;               % 'acquired' data
xhat = A' * (d(:) .* dcf(:));   % reconstruct (perform inufft along x)

subplot(131); plot(xtrue); title('true object');
subplot(132); plot(abs(xhat)); title('reconstructed');
subplot(133); plot(angle(xhat)); title('angle');
