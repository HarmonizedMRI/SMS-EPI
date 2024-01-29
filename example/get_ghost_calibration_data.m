
%% Get odd/even echo correction parameters 

% load raw data and interpolate to Cartesian grid
fn = [datdir datfile_ghostcal '.h5'];
%draw = hmriutils.epi.loadframeraw_ge(fn, etl, np, 1);   % [nFID etl np nc]
draw = hmriutils.epi.io.readframe(fn, 1);
draw = squeeze(draw(:,:, round(np/2), :));  % we only need one shot
dc = hmriutils.epi.rampsampepi2cart(draw, kxo, kxe, nx, fov(1)*100, 'spline');

% estimate linear odd/even phase offset
x = fftshift(ifft(fftshift(dc), [], 1));  % getoephase expects image space
verbose = true;
a = hmriutils.epi.getoephase(x, verbose); % 'a' contains constant and linear odd/even phase offsets

