function d2d = gedatreshape(dat, kx, npre, ntrap, ny)
%
% Inputs
% dat   [nfid ncoils]   Data for one whole EPI train
% kx    [nfid]          (cycles/cm) kspace sampling locations
%
% Usage example:
% >> [~,gx] = toppe.readmod('readout.mod'); % gx: G/cm
% >> gamma = 4.2576;      % kHz/Gauss
% >> dt = 4e-3;           % msec
% >> kx = cumsum(gx)*gamma*dt;   % cycles/cm
% >> dat = loadpfile('P_epi.7');  % [nfid ncoils nslices 1 nframes]
% >> slice = 20; frame = 30;
% >> d2d = gedatreshape(dat(:,:,slice,1,frame), kx, npre, ntrap, ny);

d2d = zeros(nro,ny);
j2d = zeros(nro,ny);
for echo = 1:ny   % EPI echo (not dabecho)
   istart = npre + (echo-1)*nro + 1;
   istop = istart + nro - 1;
   d2d(:,echo) = dat(istart:istop);
   kx2d(:,echo) = kx(istart:istop);
end
kxo = kx2d(:,1);
kxe = kx2d(:,2);

return


