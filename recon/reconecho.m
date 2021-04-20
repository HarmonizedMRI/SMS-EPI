function x = reconecho(dat, kx, nx, fov, gx)
% 1D nufft reconstruction
%
% dat: [nt 1]   
% kx:  [nt 1]  (cycles/cm)
% nx   int
% fov  cm
% gx   [nt 1]  (a.u.) readout gradient (typically a trapezoid) for one echo in EPI train

kmax = max(kx);
datg = interp1(kx, dat, linspace(kx(1),kx(end),nx));
if kx(1) > kx(end)
	datg = flipdim(datg,2);
end
datg(isnan(datg)) = 0;
x = fftshift(ifft(fftshift(datg(:),1), [], 1),1);
return;

% the nufft way:

% System matrix
xinit = zeros(nx,1);
nufft_args = {[nx],[6],[2*nx],[nx/2],'minmax:kb'};
mask = true(nx,1);
A = Gmri([fov(1)*kx(:)],mask,'nufft',nufft_args);

% Density compensation function
dcf = gx(:)/max(gx);  % simply the gradient (kspace velocity)

% recon (1d)
x = A' * (dat(:) .* dcf(:));

return;

