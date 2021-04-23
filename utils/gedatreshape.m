function [d2d, kxo, kxe] = gedatreshape(dat, kx, npre, ntrap, ny)
% function [d2d, kxo, kxe] = gedatreshape(dat, kx, npre, ntrap, ny)
%
% TOPPE acquires data continuously into one long EPI vector,
% that this function reshapes into a 2D matrix.
%
% Inputs:
% dat    [nfid ncoils nslices 1 nframes] 
% kx     [nfid]  (cycles/cm) kspace sampling locations for one EPI train
% npre   int     Number of acquired samples before first echo
% ntrap  int     Number of samples in each single-echo trapezoid
% ny     int     Number of y phase-encodes (matrix size)
%
% Outputs:
% d2d    [ntrap ny ncoils nslices nframes]  Reconstruct each coil with recon2depi.m
% kxo    [ntrap]    Odd echo kspace sampling locations
% kxe    [ntrap]    Even echo kspace sampling locations
%
% Usage example:
% >> [~,gx] = toppe.readmod('readout.mod'); % gx: G/cm
% >> gamma = 4.2576;      % kHz/Gauss
% >> dt = 4e-3;           % msec
% >> kx = cumsum(gx)*gamma*dt;   % cycles/cm
% >> dat = loadpfile('P_epi.7');  % [nfid ncoils nslices 1 nframes]
% >> d2d = gedatreshape(dat, kx, npre, ntrap, ny);

[nfid ncoils nslices nechoes nframes] = size(dat);

toppe.utils.textprogressbar('Reshaping data: ')
d2d = zeros(ntrap,ny,ncoils,nslices,nframes);
k2d = zeros(ntrap,ny);
for ic = 1:ncoils
	toppe.utils.textprogressbar(ic/ncoils*100);
	for isl = 1:nslices
		for ifr = 1:nframes
			for echo = 1:ny   % EPI echo (not dabecho)
				istart = npre + (echo-1)*ntrap + 1;
				istop = istart + ntrap - 1;
				d2d(:,echo,ic,isl,ifr) = dat(istart:istop,ic,isl,1,ifr);
				kx2d(:,echo) = kx(istart:istop);
			end
		end
	end
end
kxo = kx2d(:,1);
kxe = kx2d(:,2);
toppe.utils.textprogressbar(' done.');

return


