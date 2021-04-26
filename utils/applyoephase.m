function [d2d, kxo, kxe] = applyoephase(ph, d2d, kxo, kxe)
%
% Inputs
%  ph     [nslices 2]   odd/even phase fit parameters, see getoephase.m
%                       ph(:,1)   constant phase offset (rad)
%                       ph(:,2)   linear phase offset (rad/fov)
%  d2d    [ntrap ny ncoils nslices nframes]  See gedatreshape.m
%  kxo:   [ntrap]   (cycles/cm) sampling locations for odd echoes
%  kxe:   [ntrap]   (cycles/cm) sampling locations for even echoes
%
% Outputs:
%  d2d    same as input except after subtracting constant odd/even phase
%  kxosl  [ntrap nslices]  kxo after shifting by ph(:,2)/(2*pi) samples 
%  kxesl  [ntrap nslices]  kxe                  "

[ntrap ny ncoils nslices nframes] = size(d2d);


ntrap = length(kxo);
k.x = zeros(ntrap, ny);
k.y = zeros(ntrap, ny);

kxosl = zeros(ntrap, nslices);
kxesl = zeros(ntrap, nslices);

for isl = 1:nslices
	% apply constant phase offset
	d2d(:,1:2:end,:,isl,:) = exp(+1i*ph(isl,1)/2)*d2d(:,1:2:end,:,isl,:);
	d2d(:,2:2:end,:,isl,:) = exp(-1i*ph(isl,1)/2)*d2d(:,2:2:end,:,isl,:);

	% apply k-space shift
	ktmp = [kxo(1)*ones(4,1); kxo; kxo(end)*ones(4,1)]; % to avoid NaN after interpolation
	tmp = interp1(1:length(ktmp), ktmp, (1:length(ktmp)) + ph(isl,2)/2/(2*pi));
	kxosl(:,isl) = tmp(5:(end-4));

	ktmp = [kxe(1)*ones(4,1); kxe; kxe(end)*ones(4,1)];
	tmp = interp1(1:length(ktmp), ktmp, (1:length(ktmp)) + ph(isl,2)/2/(2*pi));
	kxesl(:,isl) = tmp(5:(end-4));
end

return
