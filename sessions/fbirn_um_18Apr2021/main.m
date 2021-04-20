
% sensitivity maps (and later: B0 field map)
%d = toppe.utils.loadpfile('~/tmp/smsepi/P_gre3d.7');
if 0
	pfile = 'P_gre3d.7';
	system('tar xf gre3d.tar readout.mod');
	[ims imsos d] = toppe.utils.recon3dft(pfile, ...   % size(d) = [256    64    64    32] 
		'readoutFile', 'readout.mod', ...
		'echo', 1);
	%tic; sens_bart = bart('ecalib', d(2:4:end,:,:,:)); toc;  % takes 19 min
	d = d(2:4:end,:,:,:);  % oprbw = 31.25 kHz (decimation = 4)
	fprintf('getting bart sens maps...');
	tic; sens_bart = bart('ecalib -r 20', d); toc;   % takes 13.6 min
	fprintf('\n');
	save sens_bart sens_bart
else
	load sens_bart
	sens = sens_bart(:,:,:,:,1);
	clear sens_bart
	sens = flipdim(sens,1);   % bart seems to flip the first dim
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
IZ = [ nz/2-isl, ny/2, nz/2+isl-1, nz-5];
sens = sens(:,:,IZ,:);

% blipped CAIPI sampling pattern
IZ = caipi(ny,3,1);  % NB! Acquisition was mb=3, so kz=4 not sampled (recon requires mb=even)

% image support
imask = true(imsize);
imask(:,:,end) = false;

% acquired data
dataprep;  % dcart

% reconstruct
xhat = reconsms(dcart(:), IZ, imask, sens, 1e-6);
im(xhat);

return
