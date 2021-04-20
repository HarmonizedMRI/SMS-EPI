% toy object, but using real sens maps

% mb factor
mb = 4;

% sensitivity maps
load sens_bart;
sens = sens_bart(:,:,6:10:(mb*10+4),:,1);
clear sens_bart
sens = flipdim(sens,1);   % bart seems to flip the first dim
ncoils = size(sens,4);

% object
imsize = [64 64 mb];
n = imsize(1);
nz = imsize(3);
fov = 20;  % cm
clear xtrue
for iz = 2:(nz-1)
	xtrue(:,:,iz) = phantom(n) * (-1)^(iz+1) * iz/nz;
end
xtrue(n/4:3*n/4,n/4:3*n/4,1) = 1;
xtrue(n/4:3*n/4,n/4:3*n/4,iz+1) = 0.5;
for iz = 1:nz
	xtrue(:,:,iz) = imrotate(xtrue(:,:,iz), 90*(iz-1));
end
%xtrue(:,:,end) = 0;

[nx ny nz] = size(xtrue);

% object support
ss = sqrt(sum(abs(sens).^2,4));
imask = ss > 0.05*max(ss(:));
imask = true(imsize);

% blipped CAIPI sampling pattern
skip = 2;
IZ = caipi(n,mb,skip);

% readout gradient (trapezoidal gradient)
dt = 4e-6;             % gradient raster time (s)
gamma = 4257.6;        % Hz/G
res = fov/nx;          % spatial resolution (cm)
kmax = 1/(2*res);      % cycles/cm
area = kmax/gamma;     % G/cm * sec (area of each readout trapezoid)
gmax = 1/(fov*gamma*dt);    % Gauss/cm
gslew = 10;      % G/cm/ms
gx = toppe.utils.trapwave2(2*area, gmax, gslew, dt*1e3);
gx = gx(2:(end-1));
kx = gamma*dt*cumsum(gx);
kx = kx - max(kx)/2;  % cycles/cm

% synthesize noisy test sms data with ramp sampling
nufft_args = {[nx],[6],[2*nx],[nx/2],'minmax:kb'};
mask = true(nx,1);
A = Gmri([fov*kx(:)],mask,'nufft',nufft_args);  % for ramp sampling
fprintf('Synthesize ramp sampled EPI data... ');
tic;
%A = getAsms(IZ, imask, sens);
%y = A*xtrue(imask);   % this should be the same as the following loop
clear d2d;
for ic = 1:ncoils
   tmp = xtrue.*sens(:,:,:,ic);
   tmp = fftshift(fft(fftshift(tmp,2), [], 2), 2);
   tmp = fftshift(fft(fftshift(tmp,3), [], 3), 3);
	for iy = 1:ny
		d2d(:,iy,ic) = A*tmp(:,iy,IZ(iy));
	end
end
toc;

SNR = 4;
d2d = d2d + randn(size(d2d))*mean(abs(d2d(:)))/SNR;

% interpolate onto cartesian grid along readout
clear dcart;
[~,A,dcf] = reconecho([], nx, [], [], kx, fov);
tic;
fprintf('Interpolate back onto cartesian grid... ');
for ic = 1:ncoils
	x = zeros(nx,ny);
	for iy = 1:ny
		x(:,iy) = reconecho(d2d(:,iy,ic), nx, A, dcf);
	end
	dcart(:,:,ic) = fftshift(fft(fftshift(x,1), [], 1),1);
end
dcart = reshape(dcart, [], ncoils);  % [sum(kmask(:)) ncoils]
toc;

% reconstruct
tol = 1e-5;
fprintf('Reconstructing...\n');
xhat = reconsms(dcart(:), IZ, imask, sens, tol);
im(xhat); colormap jet; 

return;




load tmp/info   % gx1
kx1 = cumsum(gx1);
kx1 = kx1/max(kx1(:)) - 0.5;
for ic = 1:ncoils
	tmp = fftshift(fftn(fftshift(sens(:,:,:,ic).*xtrue)));  % [nx ny nz]
	for iy = 1:ny
		tmpr = interp1(linspace(kx1(1), kx1(end), nx), tmp(:,iy,IZ(iy)), kx1);
		tmpr(isnan(tmpr)) = 0;
		if kx1(1) > kx1(end)
   		tmpr = flipdim(tmpr,2);
 		end
		if mod(iy,2)
			d2d(:,iy,ic) = flipdim(tmpr,2);
			kx2d(:,iy) = flipdim(kx1,2);
		else
			d2d(:,iy,ic) = tmpr;
			kx2d(:,iy) = kx1;
		end
			d2d(:,iy,ic) = tmpr;
			kx2d(:,iy) = kx1;
 	end
end


% old

type = 'leak';
nbrs = 4;
chat = 0;
dist_power = 1;
kappa = arg.imask;
[C, ~] = C2sparse(type, kappa, nbrs, chat, dist_power);
W = Gdiag(ones(size(A,1),1));   % weighting matrix
C = 0; %Gdiag(zeros(arg.np,1));
xinit = zeros(imsize);
tic; [xhat, info] = qpwls_pcg1(xinit(:), A, W, y, C, 'niter', 100); toc;
