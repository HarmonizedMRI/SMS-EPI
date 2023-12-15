
%% Get odd/even echo correction parameters 

% load raw data and interpolate to Cartesian grid
fn = [datdir pfile_cal]; 
draw = hmriutils.epi.loadframeraw_ge(fn, etl, np, 1);   % [nFID etl np nc]
draw = squeeze(draw(:,:, round(np/2), :));  % only need one shot
dc = hmriutils.epi.rampsampepi2cart(draw, kxo, kxe, nx, fov(1)*100, 'spline');

% estimate linear odd/even phase offset
x = fftshift(ifft(fftshift(dc), [], 1));  % getoephase expects image space
verbose = false;
a = hmriutils.epi.getoephase(x, verbose); % 'a' contains constant and linear odd/even phase offsets
save a a

%% Get single-slice EPI data for slice GRAPPA calibration

% load raw data
fn = [datdir pfile_mb1]; 
draw = hmriutils.epi.loadframeraw_ge(fn, etl, nz, 1, true);   % [nFID etl np nc]
dcal = hmriutils.epi.rampsampepi2cart(draw, kxo, kxe, nx, fov(1)*100, 'spline');

% apply odd/even phase correction
dcal = hmriutils.epi.epiphasecorrect(dcal, a);

% undo slice interleaving (see writeEPI.m)
Z = [1:2:nz 2:2:nz];
Z = [Z(1:end-2) Z(end) Z(end-1)]; % last two shots/partitions are swapped
dcal(:,:,Z,:) = dcal;

% recon and display
Icalrss = zeros(nx, etl, nz);
for iz = 1:nz
    [~, Icalrss(:,:,iz)] = toppe.utils.ift3(squeeze(dcal(:,:,iz,:)), 'type', '2d');
end
im(Icalrss);


