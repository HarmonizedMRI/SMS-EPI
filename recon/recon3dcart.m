function xhat = recon3dcart(dat, kmask, imask, sens)
% function xhat = recon3dcart(dat, kmask, imask, sens)
%
% Undersampled 3d cartesian reconstruction.
%
% Uses iterative recon for convenience and flexibility: image support, 
% and can add regularization and b0/T2* maps later
%
% Inputs:
%  dat    [nt*ncoils 1]      acquired data (complex). nt = sum(kmask(:)).
%  kmask  [nx ny nz]         logical mask indicating sampled locations
%  sens   [nx ny nz ncoils]  sensitivity maps
%  imask  [nx ny nz]         image support (logical)
%
% Output:
%  xhat   [nx ny nz]  reconstructed image

A = getAsense(kmask, imask, sens);

%W = Gdiag(ones(size(A,1),1));   % weighting matrix
%C = 0; %Gdiag(zeros(arg.np,1));
xinit = zeros(size(imask));
%tic; [xhat, info] = qpwls_pcg1(xinit(:), A, W, y, C, 'niter', 100); toc;
tol = 1e-4; nitmax = 200;
tic; [xhat,res] = cgnr_jfn(A, dat, xinit(imask), nitmax, tol); toc;
xhat = embed(xhat, imask);

return;

% compare with 
x2 = reshape(A'*y, arg.imsize);
%im(cat(1,x2,xhat)); colormap jet; 
subplot(121); im(x2); colormap jet; 
subplot(122); im(xhat); colormap jet; 

return

