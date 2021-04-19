
% sensitivity maps (and later: B0 field map)
%d = toppe.utils.loadpfile('~/tmp/smsepi/P_gre3d.7');
pfile = '~/tmp/smsepi/P_gre3d.7';
[ims imsos] = toppe.utils.recon3dft(pfile, 'readoutFile', '~/tmp/smsepi/readout_gre3d.mod', 'echo', 1);

