function x = reconecho(dat, kx, nx, fov, gx)
% 1D nufft reconstruction
%
% dat: [nt 1]   
% kx:  [nt 1]  (cycles/cm)
% nx   int
% fov  cm
% gx   [nt 1]  readout gradient (typically a trapezoid) for one echo in EPI train

% System matrix
xinit = zeros(nx,1);
nufft_args = {[nx],[6],[2*nx],[nx/2],'minmax:kb'};
mask = true(nx,1);
A = Gmri([fov(1)*kx(:)],mask,'nufft',nufft_args);

% Density compensation function
dcf = gx(:)/max(gx);  % simply the gradient (kspace velocity)

% recon (1d)
x = A' * (dat(:) .* dcf(:));   % reconstruct (perform inufft along x)

return;

