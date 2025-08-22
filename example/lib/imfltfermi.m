function [Iout, flt] = imfltfermi(Iin, fltwidth, transitionwidth)
% function [Iout, flt] = imfltfermi(Iin, fltwidth, transitionwidth)
%
% Apply low-pass 1D/2D/3D Fermi filter (rectangular) to 3D image volume(s)
%
% Iin                 [nx ny nz nframes]  Image volume time-series 
% fltwidth            [1] or [3]  Filter width, percent of full k-space extent
% transitionwidth     [1] or [3]  Fermi filter roll-off (75%-25% width), percent of full k-space extent
%
% Example usage: see sub_demo.m below. To run:
%   >> imfltfermi('demo');

import toppe.utils.*

if strcmp(lower(Iin), 'demo') | strcmp(lower(Iin), 'test')
    sub_demo;
    return;
end

assert(any(ndims(Iin) == [3 4]), 'Input image must be 3 or 4 dimensional');
assert(any(length(fltwidth) == [1 3]), 'fltwidth must be a scalar or a length-3 vector');
assert(any(length(transitionwidth) == [1 3]), 'transitionwidth must be a scalar or a length-3 vector');

nframes = size(Iin,4);
N = size(Iin(:,:,:,1));

W = round(N.*fltwidth/100);
trW = round(N.*transitionwidth/100);
flt = sub_fermi3d(N, W, trW);

Iout = 0*Iin;

for ii = 1:nframes
    d = flt .* fftshift(fftn(fftshift(Iin(:,:,:,ii))));
    Iout(:,:,:,ii) = fftshift(ifftn(fftshift(d)));
end

return

function f = sub_fermi3d(N,W,trW)
% 
%  N      [3]   matrix size 
%  W      [3]   width of filter (pixels)
%  trW    [3]   25-75% transition width
%
% Based on fermi2d, Doug Noll

assert(length(N)==3, 'Matrix size (N)  must be a vector of length 3');
assert(length(W)==3, 'Filter width (W) must be a vector of length 3');
assert(all(W<=N), 'filter width cannot exceed matrix size');
assert(length(trW)==3, 'Filter transition width (trW) must be a vector of length 3');

%[W(:) trW(:)]

cent = N/2 + 1;
x = (1:N(1));
y = (1:N(2))';
z = (1:N(3))';
[X,Y,Z] = meshgrid(x, y, z);

X = X-cent(1);
Y = Y-cent(2);
Z = Z-cent(3);

fx = 1 ./ (1 + exp( (abs(X) - W(1)/2) / max(trW(1)/2,eps) ));
fy = 1 ./ (1 + exp( (abs(Y) - W(2)/2) / max(trW(2)/2,eps) ));
fz = 1 ./ (1 + exp( (abs(Z) - W(3)/2) / max(trW(3)/2,eps) ));

if W(1) == N(1)
    fx = ones(size(fx));
end
if W(2) == N(2)
    fy = ones(size(fy));
end
if W(3) == N(3)
    fz = ones(size(fz));
end

f = fx.*fy.*fz;

return

function sub_demo

    % object
    N = [100 100 100];
    I = repmat(phantom(N(1), N(2)), [1 1 N(3)]);
    I(:,:,[1:27 end-27:end]) = 0;

    %  3D apodization, isotropic filter
    subplot(231);
    fltwidth = 80;        % percent 
    transitionwidth = 5;  % percent
    [I_filtered, flt] = imfltfermi(I, fltwidth, transitionwidth);
    im(cat(1, I(:,:,end/2), squeeze(I_filtered(:,:,end/2))));
    title('isotropic filter');

    subplot(232);
    plot(flt(:,end/2,end/2)); title(sprintf('k-space filter (width = %d%%, transitionwidth = %d%%)', fltwidth, transitionwidth));
    xlabel('pixel');
    axis([1 size(flt,1) 0 1.1]);

    subplot(233);
    x = 29; y = 50;
    plot(squeeze(I_filtered(x,y,:)), 'o-');
    title(sprintf('z profile at (x,y) = (%d,%d)', x, y));
    xlabel('pixel');

    %  anisotropic 2D apodization
    subplot(234);
    fltwidth = [50 80 100];      % percent 
    transitionwidth = [7 3.5 0];   % percent
    [I_filtered, flt] = imfltfermi(I, fltwidth, transitionwidth);
    im(cat(1, I(:,:,end/2), squeeze(I_filtered(:,:,end/2))));
    title('anisotropic 2d filter');

    subplot(235);
    title('k-space filter');
    hold off; 
    plot(squeeze(flt(:,end/2,end/2))); 
    hold on; 
    plot(squeeze(flt(end/2,:,end/2))); 
    plot(squeeze(flt(end/2,end/2,:))); 
    legend('x', 'y', 'z');
    axis([1 size(flt,1) 0 1.1]);

    subplot(236);
    plot(squeeze(I_filtered(x,y,:)), 'o-');
    title(sprintf('z profile at (x,y) = (%d,%d)', x, y));
    xlabel('pixel');

return
