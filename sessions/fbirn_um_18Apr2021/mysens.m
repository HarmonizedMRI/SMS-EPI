pfile = 'P_gre3d.7';
system('tar xf gre3d.tar readout.mod');
[ims imsos d] = toppe.utils.recon3dft(pfile, ...   % size(d) = [256 64 64 32] = [nfid ny nz ncoils]
	'readoutFile', 'readout.mod', ...
	'flipfid', true, ...
	'flipim', false, ...
	'echo', 1);

ncoils = size(ims,4);

load sens_bart
sens_bart = flipdim(sens_bart,3);
imcomb = sum(ims./sens_bart,4); % complex coil combination. Mag is very similar to imsos.

sens = ims./repmat(imsos,[1 1 1 ncoils]);

for ic = 1:4:ncoils
	figure; im(cat(1, sens(:,:,:,ic), sens_bart(:,:,:,ic)));
end
