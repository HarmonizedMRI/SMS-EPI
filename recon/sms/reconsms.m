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
arg.kmask = false(imsize);
for iz = 1:mb
	arg.kmask(:,iz:mb:end,iz) = true;
end

% 3D cartesian multi-coil system matrix
% Later: add B0 field map
%k = [kx(:) ky(:) kz(:)];
arg.imsize = imsize;
arg.imask = true(imsize);
arg.sens = sens;

arg.nc = size(arg.sens,4);
arg.nt = sum(arg.kmask(:));    % number of acquired samples (per coil)
arg.np = prod(arg.imsize);     % number of spatial positions (voxels)

if 0
	A0 = getA(arg);
	A = Asense(A0, sens);  % need to learn how to use this
else
	A = getAsense(arg);    % explicitly implements forw and back using sens maps
end

% 'acquired' data
yfull = A*xtrue(:);
SNR = 20;
yfull = yfull + randn(size(yfull))*mean(abs(yfull(:)))/SNR;

y = yfull;

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

% reconstruct
type = 'leak';
nbrs = 4;
chat = 0;
dist_power = 1;
kappa = arg.imask;
%[C, ~] = C2sparse(type, kappa, nbrs, chat, dist_power);

W = Gdiag(ones(size(A,1),1));   % weighting matrix
xinit = zeros(imsize);
%tic; [xhat, info] = qpwls_pcg1(xinit(:), A, W, y, 0, 'niter', 100); toc;
nitmax = 200;
tol = 1e-5;
tic; [xhat,res] = cgnr_jfn(A, y, xinit(:), nitmax, tol); toc;
xhat = reshape(xhat, [arg.imsize]);
im(xhat); colormap jet; 
%xcp = A'*y;
%xcp = reshape(xcp, [arg.imsize]);
%figure; im(cat(1, xcp/4, xhat)); colormap jet; 

return


