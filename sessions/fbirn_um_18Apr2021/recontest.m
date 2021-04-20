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
imask = true(imsize);
imask = ss > 0.05*max(ss(:));

% blipped CAIPI undersampling pattern
skip = 2;
IZ = caipi(n,mb,skip);

% synthesize noisy test sms data with ramp sampling
%A = getAsms(IZ, imask, sens);
%y = A*xtrue(imask);   % this should be the same as the following loop
clear y;
for ic = 1:ncoils
   tmp = fftshift(fftn(fftshift(xtrue.*sens(:,:,:,ic))));
   for iy = 1:n
      y(:,iy,ic) = tmp(:,iy,IZ(iy));
   end
end

load tmp/info   % gx1
kx1 = cumsum(gx1);
kx1 = kx1/max(kx1(:)) - 0.5;
for ic = 1:ncoils
	tmp = fftshift(fftn(fftshift(sens(:,:,:,ic).*xtrue)));  % [nx ny nz]
	for iy = 1:ny
		tmpr = interp1(linspace(kx1(1), kx1(end), nx), tmp(:,iy,IZ(iy)), kx1);
		tmpr(isnan(tmpr)) = 0;
		if kx(1) > kx(end)
   		tmpr = flipdim(tmpr,2);
 		end
		d2d(:,iy,ic) = tmpr;
		kx2d(:,iy) = kx1;
 	end
end
SNR = 4;
d2d = d2d + randn(size(d2d))*mean(abs(d2d(:)))/SNR;

% interpolate onto cartesian grid along readout
if 0
for ic = 1:ncoils
	x = zeros(nx,ny);
	for iy = 1:ny
		%x(:,iy) = reconecho(d2d(:,iy,ic), kx2d(:,iy), nx, fov, gx1); 
		x(:,iy) = interp1(kx2d(:,iy), d2d(:,iy,ic), linspace(kx2d(1,iy), kx2d(end,iy), nx)');
	end
	%dcart(:,:,ic) = fftshift(fft(fftshift(x,1), [], 1),1);
	dcart(:,:,ic) = x;
end
dcart = reshape(dcart, [], ncoils);  % [sum(kmask(:)) ncoils]
end

% reconstruct
tol = 1e-5;
xhat = reconsms(y(:), IZ, imask, sens, tol);
im(xhat); colormap jet; 

return;


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
