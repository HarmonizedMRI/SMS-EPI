% toy object, but using real sens maps 

% mb factor
mb = 6;

% define SMS slices
slSep = 2;  % slice separation (cm)
Z = [(-mb/2+0.5):(mb/2-0.5)]*slSep;  % slice locations (cm)

% sensitivity maps
load sens_bart;  % [64 64 64 32]
sens = sens_bart;  clear sens_bart
zfov = 25.6;  % cm
zres = zfov/size(sens,3);
zind = size(sens,3)/2 + round(Z/zres); % SMS slices
senstrue = sens(:,:,zind,:);  
sens = sens(:,:,zind+0,:);  % test impact of wrong sens map slices

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
xtrue = xtrue.*exp(1i*pi/2*xtrue);  % make it complex

[nx ny nz] = size(xtrue);

% object support
ss = sqrt(sum(abs(sens).^2,4));
imask = ss > 0.05*max(ss(:));
imask = true(imsize);

% synthesize noisy EPI sms data (Cartesian)
skip = 1;
IZ = caipi(n,mb,skip);
kzmax = 1/(2*slSep); % cycles/cm
KZ = (IZ-mb/2-0.5)/(mb/2)*kzmax; 
y = zeros(nx, ny, ncoils);
for ic = 1:ncoils
	for iy = 1:ny
		x = 0*xtrue;
		for iz = 1:mb
			x(:,:,iz) = exp(1i*2*pi*KZ(iy)*Z(iz)) * senstrue(:,:,iz,ic) .* xtrue(:,:,iz);
		end
		xsum = sum(x,3);
		tmp = fftshift(fftn(fftshift(xsum)));
		y(:,iy,ic) = tmp(:,iy);
	end
end
y = y + randn(size(y))*mean(abs(y(:)))/3;

% reconstruct
fprintf('Reconstructing...\n');
A = Gsms(KZ, Z, sens, imask);
%y = A*xtrue(:);
xinit = zeros(size(imask));
tol = 1e-6; nitmax = 15;
tic; [xhat,res] = cgnr_jfn(A, y(:), xinit(imask), nitmax, tol); toc;
%W = 1; C = 0;
%xhat = qpwls_pcg1(xinit(imask), A, W, y(:), C, 'niter', 250);  % runs but doesn't find the right solution
xhat = embed(xhat, imask);
im(xhat)
%subplot(121); im(xhat)
%subplot(122); plot(res, 'o-');

return;
