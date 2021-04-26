function d2d = applyoephase(d2d, ph)
% Add odd/even phase to 2D cartesian data to emulate real acquisition.
%
% Inputs
%  d2d    [nx ny nslices]  Cartesian data
%  ph     [nslices 2]   odd/even phase fit parameters, see getoephase.m
%                       ph(:,1)   constant phase offset (rad)
%                       ph(:,2)   linear phase offset (rad/fov) (gradient/acquisition delay)
%
% Outputs:
%  d2d    same as input except after adding constant and linear odd/even phase

[nx ny nslices] = size(d2d);

for isl = 1:nslices
	% apply constant phase offset
	d2d(:,1:2:end,isl) = exp(+1i*ph(isl,1)/2)*d2d(:,1:2:end,isl);
	d2d(:,2:2:end,isl) = exp(-1i*ph(isl,1)/2)*d2d(:,2:2:end,isl);

	% apply shift
	for iy = 1:2:ny
   	tmp = [d2d(1,iy)*ones(4,1); d2d(:,iy); d2d(end,iy)*ones(4,1)]; % to avoid NaN after interpolation
   	tmp = interp1(1:size(tmp), tmp, (1:length(tmp)) - ph(isl,2)/2/(2*pi));
   	d2d(:,iy) = tmp(5:(end-4));
	end
	for iy = 2:2:ny
   	tmp = [d2d(1,iy)*ones(4,1); d2d(:,iy); d2d(end,iy)*ones(4,1)]; % to avoid NaN after interpolation
   	tmp = interp1(1:size(tmp), tmp, (1:length(tmp)) + ph(isl,2)/2/(2*pi));
   	d2d(:,iy) = tmp(5:(end-4));
	end
end

return
