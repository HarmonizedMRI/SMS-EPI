
% sensitivity maps (and later: B0 field map)
%d = toppe.utils.loadpfile('~/tmp/smsepi/P_gre3d.7');
if 0
pfile = 'P_gre3d.7';
system('tar xf gre3d.tar readout.mod');
[ims imsos d] = toppe.utils.recon3dft(pfile, ...   % size(d) = [256    64    64    32] 
	'readoutFile', 'readout.mod', ...
	'echo', 1);
%tic; sens_bart = bart('ecalib', d(2:4:end,:,:,:)); toc;  % takes 19 min
d = d(2:4:end,:,:,:);  % oprbw = 31.25 kHz (decimation = 4)
fprintf('getting bart sens maps...');
tic; sens_bart = bart('ecalib -r 20', d); toc;   % takes 13.6 min
fprintf('\n');
sens_bart = sens_bart(:,:,:,:,1);
save sens_bart sens_bart
else
load sens_bart
sens = sens_bart(:,:,:,:,1);
clear sens_bart
sens = flipdim(sens,1);   % bart seems to flip the first dim
end

% matrix size for reconstruction
mb = 3; % multiband/sms factor (number of simultaneous slices)
[nx ny] = size(sens(:,:,1,1));
imsize = [nx ny mb];

% blipped CAIPI sampling pattern
kmask = false(imsize);
for iz = 1:mb
	kmask(:,iz:mb:end,iz) = true;
end

% image support
imask = true(imsize);

% acquired data
pfile = 'P_smsepi.7';
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile); % dat = [8292 32 1 1 40] (fid ncoils nslices nechoes nviews)
dat = flipdim(dat,1); % yikes
frame = 10;
dat = dat(:,:,1,1,frame);  % [nfid ncoils]
ncoils = size(dat,2);

% get scan info
load tmp/info  % gpre gx1 kx fov (run smsepi.m)

% EPI correction parameters
delay = 0.16;  % fraction of 4us sample
th0 = 0.2;   % odd/even dc phase offset

% apply temporal shift (odd/even linear phase correction)
nt = length(kx);
kx = interp1(1:nt, kx, (1:nt)-delay);

% reshape 
for echo = 1:ny   % EPI echo (not dabecho)
	istart = length(gpre) + (echo-1)*length(gx1) + 1;
	istop = istart + length(gx1) - 1;
	d2d(:,:,echo) = dat(istart:istop, :);
	kx2d(:,echo) = kx(istart:istop);
end
d2d = permute(d2d, [1 3 2]);  % [length(gx1) 64 ncoils]

% apply odd/even dc phase offset
d2d(:,2:2:end) = bsxfun(@times, exp(1i*th0), d2d(:,2:2:end));

% interpolate onto cartesian grid along readout
for ic = 1:ncoils
	fprintf('%d\n', ic);
	for iy = 1:ny
		x(:,iy) = reconecho(d2d(:,iy,ic), kx2d(:,iy), nx, fov, gx1); 
	end
	dcart(:,:,ic) = fftshift(fft(fftshift(x,1), [], 1),1);
end
dcart = reshape(dcart, [], ncoils);  % [sum(kmask(:)) ncoils]

% reconstruct
%xhat = recon3dcart(dcart, kmask, imask, sens);
%im(xhat); colormap jet; 

return
