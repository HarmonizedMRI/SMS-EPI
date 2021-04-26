function d2d = applyoephase(d2d, ph)
% Add odd/even phase to 2D cartesian data to emulate real acquisition.
%
% Inputs
%  d2d   [nx ny nslices]  Cartesian data
%  ph    [nslices 3]   odd/even phase fit parameters, see getoephase.m
%                      ph(:,1)   constant phase offset (rad)
%                      ph(:,2)   linear phase offset along x (rad/fov) (gradient/acquisition delay)
%                      ph(:,3)   linear phase offset along y
%
% Outputs:
%  d2d    same as input except after adding constant and linear odd/even phase

[nx ny nslices] = size(d2d);

for isl = 1:nslices
	% apply constant phase offset
	d2d(:,1:2:end,isl) = exp(-1i*ph(isl,1)/2)*d2d(:,1:2:end,isl);
	d2d(:,2:2:end,isl) = exp(+1i*ph(isl,1)/2)*d2d(:,2:2:end,isl);

	% apply shift in x 
	for iy = 1:1:ny
   	tmp = [d2d(1,iy)*ones(4,1); d2d(:,iy); d2d(end,iy)*ones(4,1)]; % to avoid NaN after interpolation
   	tmp = interp1(1:length(tmp), tmp, (1:length(tmp)) + (-1)^(iy-1)*ph(isl,2)/2/(2*pi));
   	d2d(:,iy) = tmp(5:(end-4));
	end

	% apply shift in y
	y = 1:(nx+8);
	y(1:2:end) = y(1:2:end) + ph(isl,3)/2/(2*pi);
	y(2:2:end) = y(2:2:end) - ph(isl,3)/2/(2*pi);
	for ix = 1:nx
   	tmp = [d2d(ix,1)*ones(1,4) d2d(ix,:) d2d(ix,end)*ones(1,4)];
   	tmp = interp1(1:length(tmp), tmp, y);
   	d2d(ix,:) = tmp(5:(end-4));
	end
end

return
