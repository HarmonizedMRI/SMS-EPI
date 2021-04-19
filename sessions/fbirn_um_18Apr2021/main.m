
% sensitivity maps (and later: B0 field map)
%d = toppe.utils.loadpfile('~/tmp/smsepi/P_gre3d.7');
pfile = '~/tmp/smsepi/P_gre3d.7';
[ims imsos d] = toppe.utils.recon3dft(pfile, ...   % size(d) = [256    64    64    32] 
	'readoutFile', '~/tmp/smsepi/readout_gre3d.mod', ...
	'echo', 1);
%tic; sens_bart = bart('ecalib', d(2:4:end,:,:,:)); toc;  % takes 19 min
d = d(2:4:end,:,:,:);  % oprbw = 31.25 kHz (decimation = 4)
calsize = 20;
n = size(d,2); % matrix size (assumed isotropic)
r = (n/2-calsize/2):(n/2+calsize/2-1);
fprintf('getting bart sens maps...');
%tic; sens_bart = bart('ecalib -r 20', d); toc;   % takes 
fprintf('\n');
% save sens_bart sens_bart

