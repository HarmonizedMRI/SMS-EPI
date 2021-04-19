% toy example

% object
imsize = [64 64 6];
fov = [24 24 4];   % cm
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

% sensitivity maps
% 4 rings of 8 coils (like the Nova 32ch head coil)
nRings = 4;
nCoilsPerRing = 8;
nc = nRings*nCoilsPerRing;
[x y z] = meshgrid(linspace(-1,1,imsize(1)), linspace(-1,1,imsize(2)), linspace(-1,1,imsize(3))); 
for ir = 1:nRings
	for ic = 1:nCoilsPerRing
		p0 = sqrt(2)*exp(1i*2*pi*ic/nCoilsPerRing);
		xc = real(p0); yc = imag(p0); zc = (-nRings/2 + ir - 0.5);  % coil location
		r = sqrt((x-xc).^2 + (y-yc).^2 + (z-zc).^2);
		sens(:,:,:,(ir-1)*nCoilsPerRing+ic) = exp(-r.^2 + 1i*pi*(r-0.5)); % add some phase too
	end
end
if 0
for ic = 1:nc
	figure; subplot(121); im(sens(:,:,:,ic), [0 1]);
	subplot(122); im(angle(sens(:,:,:,ic)), [0 1]);
	colormap jet;
end
return;
end

% blipped CAIPI undersampling pattern
mb = size(xtrue,3);
kmask = false(imsize);
for iz = 1:mb
	kmask(:,iz:mb:end,iz) = true;
end

% image support
imask = true(imsize);
imask = xtrue ~= 0; %true(imsize);

% synthesize 'acquired' undersampled multicoil data
A = getAsense(kmask, imask, sens);
y = A*xtrue(imask);
SNR = 4;
y = y + randn(size(y))*mean(abs(y(:)))/SNR;

% reconstruct
xhat = recon3dcart(y, kmask, imask, sens);
im(xhat); colormap jet; 

return;




% old

% reconstruct
type = 'leak';
nbrs = 4;
chat = 0;
dist_power = 1;
kappa = arg.imask;
%[C, ~] = C2sparse(type, kappa, nbrs, chat, dist_power);

W = Gdiag(ones(size(A,1),1));   % weighting matrix
C = 0; %Gdiag(zeros(arg.np,1));
xinit = zeros(imsize);
%tic; [xhat, info] = qpwls_pcg1(xinit(:), A, W, y, C, 'niter', 100); toc;
tol = 1e-4; nitmax = 200;
tic; [xhat,res] = cgnr_jfn(A, y, xinit(:), nitmax, tol); toc;
xhat = reshape(xhat, [arg.imsize]);
x2 = reshape(A'*y, imsize);
%im(cat(1,x2,xhat)); colormap jet; 
subplot(121); im(x2); colormap jet; 
subplot(122); im(xhat); colormap jet; 

return


% misc old code

if 0   % check coil data
for ic = 1:nc
	tmp = y(((ic-1)*arg.nt+1):(ic*arg.nt));
	tmp = embed(tmp, arg.kmask);
	x = fftshift(ifftn(fftshift(tmp)));
	figure; im(cat(1, xtrue.*sens(:,:,:,ic), x))
end
return;
end

if 0
% check A_back 
xhat = A'*y(:);
xhat = reshape(xhat, [arg.imsize]);
im(cat(1, xtrue(:,:,:,1), xhat/3)); colormap jet;
return
end
