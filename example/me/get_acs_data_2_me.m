%% for Multi-echo EPI 
%% modified by Yongli He (yonglihe@umich.edu)
%% Get single-slice EPI data for slice GRAPPA calibration

% Updated 22-Nov-2024 to use frame 2 

% load raw data, interpolate to Cartesian grid, and apply odd/even phase correction
draw = hmriutils.epi.io.readframe(h5file_mb1, 2);

iecho = 2; %which echo to use for reference
draw = draw(:,etl/nTE*(iecho-1)+1:etl/nTE*iecho,:,:);

%load('kxoe80_2d.mat')
load(E.readout_trajectory_file);
nfid = numel(kxo);
kspace_delay = -3.4;
kxo = interp1(1:nfid, kxo, (1:nfid) + kspace_delay,'linear','extrap');
kxe = interp1(1:nfid, kxe, (1:nfid) + kspace_delay, 'linear','extrap'); % apply ad-hoc constant delay to avoid phase wrap in the following estimation of a

dcal = hmriutils.epi.rampsampepi2cart(draw, kxo, kxe, nx, fov(1)*100, 'spline');

%%
if 1
    deltaa = [-0.2412,0]';
    dcal = hmriutils.epi.epiphasecorrect(dcal, deltaa);
    % undo slice interleaving (see writeEPI.m)
    Z = hmriutils.epi.getsliceordering(nz);
    dcal(:,:,Z,:) = dcal;
else
    dcal = hmriutils.epi.epiphasecorrect(dcal, a);
    % undo slice interleaving (see writeEPI.m)
    Z = hmriutils.epi.getsliceordering(nz);
    dcal(:,:,Z,:) = dcal;
end

%figure;im(squeeze(abs(dcal(:,:,22,1))).^0.3,sprintf('kdelay = %.1f',kspace_delay))

% recon and display
Icalrss = zeros(nx, etl/nTE, nz);
for iz = 1:nz
    [~, Icalrss(:,:,iz)] = toppe.utils.ift3(squeeze(dcal(:,:,iz,:)), 'type', '2d');
end
figure;
im(Icalrss);

