
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
sens_bart = sens_bart(:,:,:,:,1);
save sens_bart sens_bart
else
load sens_bart
sens = sens_bart(:,:,:,:,1);
clear sens_bart
sens = flipdim(sens,1);   % bart seems to flip the first dim
end

mbfactor = 3;
[nx ny] = size(sens(:,:,1,1));
imsize = [nx ny mbfactor];

% blipped CAIPI undersampling pattern
kmask = 

% image support
imask = true(imsize);

% acquired data
% d = [nt*ncoil], where nt = sum(kmask(:))
% reconstruct
xhat = recon3dcart(y, kmask, imask, sens);
im(xhat); colormap jet; 

return
