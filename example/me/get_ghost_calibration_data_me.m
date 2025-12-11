%% Modified for multi-echo EPI 
% by Yongli He (yonglihe@umich.edu)
%% Get odd/even echo correction parameters 

% load raw data and interpolate to Cartesian grid
%draw = hmriutils.epi.loadframeraw_ge(fn, etl, np, 1);   % [nFID etl np nc]
draw = hmriutils.epi.io.readframe(h5file_ghostcal, 1);
iecho = 2; 
draw = squeeze(draw(:,(iecho-1)*etl+1:iecho*etl, round(np/2), :));  % we only need one shot, one echo 

% load('kxoe80_cal.mat')%load(E.readout_trajectory_file);
% nfid = numel(kxo);
% kspace_delay = -0.4;%-1.0;
% kxo = interp1(1:nfid, kxo, (1:nfid) + kspace_delay,'linear','extrap');
% kxe = interp1(1:nfid, kxe, (1:nfid) + kspace_delay, 'linear','extrap'); % apply ad-hoc constant delay to avoid phase wrap in the following estimation of a

dc = hmriutils.epi.rampsampepi2cart(draw, kxo, kxe, nx, fov(1)*100, 'spline');

% estimate linear odd/even phase offset
x = fftshift(ifft(fftshift(dc), [], 1));  % getoephase expects image space
verbose = true;
a = hmriutils.epi.getoephase(x, verbose); % 'a' contains constant and linear odd/even phase offsets

save a a