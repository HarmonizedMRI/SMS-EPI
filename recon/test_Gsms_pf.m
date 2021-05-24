% create toy object, simulate an SMS EPI acquisition,
% and reconstruct 

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
	xtrue(:,:,iz) = phantom(n) * iz/nz; % * (-1)^(iz+1) * iz/nz;
end
xtrue(n/4:3*n/4,n/4:3*n/4,1) = 0.25;
xtrue(n/4:3*n/4,n/4:3*n/4,iz+1) = 0.5;
for iz = 1:nz
	xtrue(:,:,iz) = imrotate(xtrue(:,:,iz), 90*(iz-1));
end

% add some phase
[xg,yg,zg] = ndgrid((0.5:n)/(n/2)-1, (0.5:n)/(n/2)-1, (0.5:nz)/(nz/2)-1);
th = exp(1i*pi/2*(xg.^2 + yg + zg.^1)); % + 1i*pi/2*0.3*abs(xtrue)/max(abs(xtrue(:))) ); 
xtrue = xtrue.*th;

[nx ny nz] = size(xtrue);

% object support
ss = sqrt(sum(abs(sens).^2,4));
imask = true(imsize);
imask = ss > 0.05*max(ss(:));

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

[nx ny nz] = size(xtrue);

% PF sampling
yfull = y;
KZfull = KZ;
y = y(:, (end/2-ny/4+1):end, :);
KZ = KZ((end/2-ny/4+1):end);

% low-res image (central fully sampled) for phase estimation. Assumes 3/4 PF.
imlo = 0*xtrue;                 % low-res image phase estimate
for iz = 1:nz
	imlo(:,:,iz) = toppe.utils.imfltfermi(xtrue(:,:,iz), nx/2, nx/4, 'circ');
%	d = fftshift(fft2(fftshift(xtrue(:,:,iz))));
%	d([1:(end/2-nx/4) (end/2+nx/4+1):end], [1:(end/2-ny/4) (end/2+ny/4+1):end]) = 0; 
%	imlo(:,:,iz) = fftshift(ifft2(fftshift(d))); 
end
imlo = exp(1i*angle(imlo));

% reconstruct
% first recon using zero-filling to get a decent initialization
fprintf('Reconstructing...\n');

A = Gsms(KZfull, Z, sens, imask); nitmax = 10;
xinit = zeros(size(imask));
tol = 1e-6;
tic; [xhat1,res1] = cgnr_jfn(A, yfull(:), xinit(imask), nitmax, tol); toc;
xhat1 = embed(xhat1, imask);

%xinit = zeros(size(imask));
xinit = xhat1 .* conj(imlo);
A = Gsms_pf(KZfull, Z, sens, imask, imlo, 'pf', 1); nitmax = 10;
tol = 1e-6;
tic; [xhat2,res2] = cgnr_jfn(A, yfull(:), xinit(imask), nitmax, tol); toc;
xhat2 = embed(xhat2, imask);

subplot(131); im(xhat1)
subplot(132); im(xhat2)
subplot(133); plot([res1 res2], 'o-');

return;
