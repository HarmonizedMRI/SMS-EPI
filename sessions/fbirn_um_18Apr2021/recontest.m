% toy example

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

% object support
ss = sqrt(sum(abs(sens).^2,4));
imask = ss > 0.05*max(ss(:));
imask = true(imsize);

% blipped CAIPI undersampling pattern
kmask = false(imsize);
for iz = 1:mb
	kmask(:,iz:mb:end,iz) = true;
end
IZ = caipi(n,mb);

% synthesize 'acquired' undersampled multicoil data
%A = getAsms(IZ, imask, sens);
%y = A*xtrue(imask);
clear y;
for ic = 1:ncoils
	tmp = fftshift(fftn(fftshift(xtrue.*sens(:,:,:,ic))));
	for iy = 1:n
		y(:,iy,ic) = tmp(:,iy,IZ(iy));
	end
end

% add noise
SNR = 4;
y = y + randn(size(y))*mean(abs(y(:)))/SNR;

% reconstruct
xhat = reconsms(y(:), IZ, imask, sens);
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
