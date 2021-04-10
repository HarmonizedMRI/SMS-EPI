% curdir = pwd; cd /opt/matlab/toolbox/irt/; setup; cd(curdir);

% run ../sequence/fmri2depi;
[kx,~] = toppe.utils.g2k([gx1(:) gx1(:)]);

% 1d recon of one echo using nufft
%npix = nx;
%xinit = zeros(npix^2,1);
xinit = zeros(nx,1);
%nufft_args = {[ny,nx],[6,6],[2*ny,2*nx],[ny/2,nx/2],'table',2^12,'minmax:kb'};
nufft_args = {[nx],[6],[2*nx],[nx/2],'minmax:kb'};
mask = true(nx); % Mask for support (image space)
L = 6;
A = Gmri([fov(1)*kx(:) fov(2)*ky(:)],mask,'nufft',nufft_args);

