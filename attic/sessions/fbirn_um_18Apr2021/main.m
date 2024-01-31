% sensitivity maps (and later: B0 field map)
if 0
	pfile = 'P_gre3d.7';
	system('tar xf gre3d.tar readout.mod');
	[ims imsos d] = toppe.utils.recon3dft(pfile, ...   % size(d) = [256 64 64 32] = [nfid ny nz ncoils]
		'readoutFile', 'readout.mod', ...
		'flipfid', true, ...
		'flipim', false, ...
		'echo', 1);

	% fix for exp on 4/18/21 (should have negated gy and gz blips to match smsepi.m)
	d = flipdim(d,2);
	d = flipdim(d,3);

	d = d(2:4:end,:,:,:);  % oprbw = 31.25 kHz (decimation = 4)
	fprintf('getting bart sens maps...');
	tic; sens_bart = bart('ecalib -r 20', d); toc;   % takes 14 min
	fprintf('\n');
	sens_bart = sens_bart(:,:,:,:,1);
	save sens_bart sens_bart
else
	load sens_bart
	sens = sens_bart;
	clear sens_bart
end

ncoils = size(sens,4);

% matrix size for reconstruction
mb = 4; % multiband/sms factor (number of simultaneous slices)
[nx ny] = size(sens(:,:,1,1));
imsize = [nx ny mb];

% pick out slices from sensitivity map
slSep = 5;  % cm (see smsepi.m)
slThick = 0.2812;
isl = round(slSep/slThick);
nz = size(sens,3);
IZmb = [ nz/2-isl, nz/2, nz/2+isl-1, nz];
sens = sens(:,:,IZmb,:);

% blipped CAIPI sampling pattern
IZ = caipi(ny,3,1);  % NB! Acquisition was mb=3, so kz=4 not sampled (recon requires mb=even)

% image support
imask = true(imsize);
imask(:,:,end) = false;

% acquired data
dataprep;  % dcart = [nx*ny ncoils]

% reconstruct
tol = 1e-6;
xhat = reconsms(dcart(:), IZ, imask, sens, tol);
im(xhat);

return
