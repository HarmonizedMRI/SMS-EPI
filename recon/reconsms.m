function xhat = reconsms(dat, IZ, imask, sens, tol, nitmax)
% function xhat = recon3dcart(dat, IZ, imask, sens)
%
% Undersampled 3d cartesian reconstruction.
%
% Uses iterative recon for convenience and flexibility: image support, 
% and can add regularization and b0/T2* maps later
%
% Inputs:
%  dat    [nt*ncoils 1]      acquired data (complex). nt = sum(kmask(:)).
%  sens   [nx ny nz ncoils]  sensitivity maps
%  imask  [nx ny nz]         image support (logical)
%
% Output:
%  xhat   [nx ny nz]  reconstructed image

if ~exist('tol', 'var')
	tol = 1e-6; 
end
if ~exist('nitmax', 'var')
	nitmax = 200;
end

A = getAsms(IZ, imask, sens);

%W = Gdiag(ones(size(A,1),1));   % weighting matrix
%C = 0; %Gdiag(zeros(arg.np,1));
xinit = zeros(size(imask));
%tic; [xhat, info] = qpwls_pcg1(xinit(:), A, W, y, C, 'niter', 100); toc;
tic; [xhat,res] = cgnr_jfn(A, dat, xinit(imask), nitmax, tol); toc;
xhat = embed(xhat, imask);

return;

% compare with 
x2 = reshape(A'*y, arg.imsize);
%im(cat(1,x2,xhat)); colormap jet; 
subplot(121); im(x2); colormap jet; 
subplot(122); im(xhat); colormap jet; 

return

