
%% Get single-slice EPI data for GRAPPA calibration


% load raw data, interpolate to Cartesian grid, and apply odd/even phase correction
draw = hmriutils.epi.io.readframe(h5file_grappacal, 1);

% iecho = 2; %which echo to use for reference
% draw = draw(:,Ry*etl*(iecho-1)+1:Ry*etl*iecho,:,:);

dcalgrappa = hmriutils.epi.rampsampepi2cart(draw, kxo, kxe, nx, fov(1)*100, 'spline');

dcalgrappa = hmriutils.epi.epiphasecorrect(dcalgrappa, a);
% undo slice interleaving (see writeEPI.m)
Z = hmriutils.epi.getsliceordering(nz);
dcalgrappa(:,:,Z,:) = dcalgrappa;



% recon and display
Ical_grappa_rss = zeros(nx, ny, nz);
for iz = 1:nz
    [~, Ical_grappa_rss(:,:,iz)] = toppe.utils.ift3(squeeze(dcalgrappa(:,:,iz,:)), 'type', '2d');
end
figure;
im(Ical_grappa_rss);

