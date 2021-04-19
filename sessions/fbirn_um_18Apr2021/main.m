
% sensitivity maps (and later: B0 field map)
%d = toppe.utils.loadpfile('~/tmp/smsepi/P_gre3d.7');
pfile = '~/tmp/smsepi/P_gre3d.7';
[ims imsos d] = toppe.utils.recon3dft(pfile, 'readoutFile', '~/tmp/smsepi/readout_gre3d.mod', 'echo', 1);
% size(d) = [256    64    64    32] 
%tic; sens2 = bart('ecalib', d(2:4:end,:,:,:)); toc;  % takes 

