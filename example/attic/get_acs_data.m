
%% Get single-slice EPI data for slice GRAPPA calibration

% load raw data, interpolate to Cartesian grid, and apply odd/even phase correction
fn = [datdir datfile_mb1 '.h5']
draw = hmriutils.epi.io.readframe(fn, 1);
dcal = hmriutils.epi.rampsampepi2cart(draw, kxo, kxe, nx, fov(1)*100, 'nufft');
dcal = hmriutils.epi.epiphasecorrect(dcal, a);

% undo slice interleaving (see writeEPI.m)
Z = hmriutils.epi.getsliceordering(nz);
dcal(:,:,Z,:) = dcal;

% recon and display
Icalrss = zeros(nx, etl, nz);
for iz = 1:nz
    [~, Icalrss(:,:,iz)] = toppe.utils.ift3(squeeze(dcal(:,:,iz,:)), 'type', '2d');
end
im(Icalrss);

