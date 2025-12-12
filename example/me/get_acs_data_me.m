%% for Multi-echo EPI 
%% modified by Yongli He (yonglihe@umich.edu)
%% Get single-slice EPI data for slice GRAPPA calibration


% load raw data, interpolate to Cartesian grid, and apply odd/even phase correction
draw = hmriutils.epi.io.readframe(h5file_mb1, 1);

iecho = 2; %which echo to use for reference
draw = draw(:,etl*(iecho-1)+1:etl*iecho,:,:);

dcal = hmriutils.epi.rampsampepi2cart(draw, kxo, kxe, nx, fov(1)*100, 'spline');

dcal = hmriutils.epi.epiphasecorrect(dcal, a);
% undo slice interleaving (see writeEPI.m)
Z = hmriutils.epi.getsliceordering(nz);
dcal(:,:,Z,:) = dcal;



% recon and display
Icalrss = zeros(nx, etl, nz);
for iz = 1:nz
    [~, Icalrss(:,:,iz)] = toppe.utils.ift3(squeeze(dcal(:,:,iz,:)), 'type', '2d');
end
figure;
im(Icalrss);

