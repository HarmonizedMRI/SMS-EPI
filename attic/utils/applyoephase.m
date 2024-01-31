function d2d = applyoephase(d2d, ph) %, b0eddy)
% function d2d = applyoephase(d2d, ph) %, b0eddy)
%
% Add odd/even phase to 2D cartesian data to emulate real acquisition.
%
% Inputs
%  d2d   [nx ny]       Cartesian data
%  ph    [nslices 3]   odd/even phase fit parameters, see getoephase.m
%                      ph(:,1)   constant phase offset (rad)
%                      ph(:,2)   linear phase offset along x (rad/fov) (gradient/acquisition delay)
%                      ph(:,3)   linear phase offset along y
%
% Outputs:
%  d2d    same as input except after adding constant and linear odd/even phase

[nx ny] = size(d2d);

npad = 4;

% apply constant phase offset
d2d(:,1:2:end) = exp(-1i*ph(1)/2)*d2d(:,1:2:end);
d2d(:,2:2:end) = exp(+1i*ph(1)/2)*d2d(:,2:2:end);

% test: add b0 eddy
%t = 1:size(d2d,1);
%t = t-t(end)/2+1;
%d2d(:,1:2:end,isl) = repmat(exp(1i*t'/max(t)*b0eddy),[1 length(1:2:nx)]).*d2d(:,1:2:end,isl);
%d2d(:,2:2:end,isl) = repmat(exp(-1i*t'/max(t)*b0eddy),[1 length(2:2:nx)]).*d2d(:,2:2:end,isl);

% apply shift in x 
for iy = 1:1:ny
	% pad to avoid nan after interpolation
  	tmp = [d2d(1,iy)*ones(npad,1); d2d(:,iy); d2d(end,iy)*ones(npad,1)]; 
  	tmp = interp1(1:length(tmp), tmp, (1:length(tmp)) + (-1)^(iy-1)*ph(2)/2/(2*pi));
  	d2d(:,iy) = tmp((npad+1):(end-npad));
end

% apply shift in y
if 1
	y = 1:(nx+2*npad);
	y(1:2:end) = y(1:2:end) + ph(3)/2/(2*pi);
	y(2:2:end) = y(2:2:end) - ph(3)/2/(2*pi);
	for ix = 1:nx
   	tmp = [d2d(ix,1)*ones(1,npad) d2d(ix,:) d2d(ix,end)*ones(1,npad)];
   	tmp = interp1(1:length(tmp), tmp, y);
   	d2d(ix,:) = tmp((npad+1):(end-npad));
	end
end

return
