function [d2d, kxo, ke] = gedatreshape(dat, kx, npre, ntrap, ny)
% function [d2d, kxo, ke] = gedatreshape(dat, kx, npre, ntrap, ny)
%
% TOPPE acquires data continuously into one long EPI vector,
% that this function reshapes into a 2D matrix.
%
% Inputs:
% dat    [nfid ncoils]   Data for one whole EPI train
% kx     [nfid]          (cycles/cm) kspace sampling locations
% npre   int             Number of acquired samples before first echo
% ntrap  int             Number of sample in each single-echo trapezoid
% ny     int             Number of y phase-encodes (matrix size)
%
% Outputs:
% d2d    [ntrap ny ncoils]   Reconstruct each coil with recon2depi.m
% kxo    [ntrap]             Odd echo kspace sampling locations
% kxe    [ntrap]             Even echo kspace sampling locations
%
% Usage example:
% >> [~,gx] = toppe.readmod('readout.mod'); % gx: G/cm
% >> gamma = 4.2576;      % kHz/Gauss
% >> dt = 4e-3;           % msec
% >> kx = cumsum(gx)*gamma*dt;   % cycles/cm
% >> dat = loadpfile('P_epi.7');  % [nfid ncoils nslices 1 nframes]
% >> slice = 20; frame = 30;
% >> d2d = gedatreshape(dat(:,:,slice,1,frame), kx, npre, ntrap, ny);

d2d = zeros(ntrap,ny);
k2d = zeros(ntrap,ny);
for echo = 1:ny   % EPI echo (not dabecho)
   istart = npre + (echo-1)*ntrap + 1;
   istop = istart + ntrap - 1;
   d2d(:,echo) = dat(istart:istop);
   kx2d(:,echo) = kx(istart:istop);
end
kxo = kx2d(:,1);
kxe = kx2d(:,2);

return


