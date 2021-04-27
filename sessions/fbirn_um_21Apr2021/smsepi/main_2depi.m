% sens maps
load ../gre3d/sens_bart
sens = sens_bart;
%sens = flip(sens,2);


%% Recon ../epi/ epi data using 2d sense nufft (for testing)

% pfile and frame to reconstruct
pfile = '../epi/P_fmri2depi.7';
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [nt ncoils nslices 1 nframes]
dat = flipdim(dat,1); % as usual

[nfid ncoils nz nframes] = size(squeeze(dat));

% kspace locations
[~,gx,gy,gz] = toppe.readmod('../epi/tmp/readout.mod'); 
gamma = 4.2576;  % kHz/G
dt = 4e-3;       % ms
k = dt*gamma*cumsum([gx gy gz]);

npre = 83;
ntrap = 101;
nx = 64; ny = 64; fov = 25.6;

% select slice and frame
frame = 8;
slice = 44;
sens = squeeze(sens(:,:,slice,:));
dat = dat(:,:,slice,1,frame);

% indeces corresponding to EPI train
istart = npre + 1;
istop = istart + ntrap*ny - 1;

% apply gradient delay
ph = [-0.14 -4.8];  % constant (rad) and linear (rad/cycle) odd/even phase offset
k(:,1) = interp1(1:nfid, k(:,1), (1:nfid) + ph(2)/(2*pi));

% crop out dephaser gradients before start of EPI train
dat = dat(istart:istop,:);
k = k(istart:istop,:);

% only use some of the coils (for testing)
%dat = dat(:,1:10:end);
%sens = sens(:,:,1:10:end);

% apply odd/even phase offset (constant)
d2d = reshape(dat, ntrap, ny, []);
for iy = 1:2:ny
  % d2d(:,1:2:end,:,isl,:) = exp(+1i*ph(isl,1)/2)*d2d(:,1:2:end,:,isl,:);
  % d2d(:,2:2:end,:,isl,:) = exp(-1i*ph(isl,1)/2)*d2d(:,2:2:end,:,isl,:);
end

% odd echoes only
ysamp = 1:2:ny;
dat = d2d(:,ysamp,:);
k = reshape(k, [ntrap ny 3]);
k = k(:,ysamp,:);
k = reshape(k, [], 3);

nufft_args = {[nx,ny],[6,6],[2*nx,2*ny],[nx/2,ny/2],'minmax:kb'};
mask = true(nx,ny); % Mask for support
L = 6;
A0 = Gmri([fov*k(:,1) fov*k(:,2)], ...
	mask, 'nufft', nufft_args);
A = Asense(A0, sens);

x0 = reshape(A'*dat(:)/(nx*ny), [nx ny]);
W = 1; C = 0;
x = qpwls_pcg1(x0, A, W, dat(:), C, ...
                   'niter', 30);
%tic; [x,res] = cgnr_jfn(A, dat(:), x0(:), 15, 1e-7); toc; % Also works
x = reshape(x, [nx ny]);
im(cat(1,x0,x));

return;


