% make sensitivity maps (and later: B0 field map)
pfile = 'P_gre3d.7';
system('tar xf gre3d.tar readout.mod');
[ims imsos d] = toppe.utils.recon3dft(pfile, ...   % size(d) = [256 64 64 32] = [nfid ny nz ncoils]
	'readoutFile', 'readout.mod', ...
	'flipfid', true, ...
	'flipim', false, ...
	'echo', 1);

d = d(2:4:end,:,:,:);  % oprbw = 31.25 kHz (decimation = 4)
fprintf('getting bart sens maps...');
tic; sens_bart = bart('ecalib -r 20', d); toc;   % takes 14 min
fprintf('\n');
sens_bart = sens_bart(:,:,:,:,1);
save sens_bart sens_bart
