function I = dcm2nii4d(nt, ofname)
% function I = dcm2nii4d(nt)
%
% Convert folder with DICOM files names 1.dcm, 2.dcm, ..., nt.dcm
% to a 4D NIFTI file.
% This is the DICOM output on Siemens in the 'Enhanced' mode (I think).
%
% nt         number of time points (3D frames)
% ofname     NIFTI output file name

d = dir('*.nii');
assert(length(d) == 0, 'Remove existing .nii files before running this command');

% Construct NIFTI header
system('dcm2niix -s y 1.dcm');
d = dir('*.nii');
niiinfo = niftiinfo(d(end).name);
[nx, ny, nz] = deal(niiinfo.ImageSize(1), niiinfo.ImageSize(2), niiinfo.ImageSize(3)); 
niiinfo.raw.dim = [4 nx ny nz nt 1 1 1];
niiinfo.ImageSize = [nx ny nz nt];
niiinfo.PixelDimensions = [niiinfo.PixelDimensions 0.8];

% Get images
I = zeros(nx, ny, nz, nt);
toppe.utils.textprogressbar('Loading images ');
for n = 1:nt
    tmp = dicomread(sprintf('%d.dcm', n));
    I(:,:,:,n) = squeeze(tmp);
    toppe.utils.textprogressbar(n/nt*100);
end
toppe.utils.textprogressbar(' ');

% Reorient so that the orientation of 'ofname' when viewed in MRIcroGL
% matches that of the .nii file created by 'dcm2niix -s y 1.dcm'
I = flipdim(permute(I, [2 1 3 4]), 2);

% Write 4D volume to .nii file
niftiwrite(int16(I), ofname, niiinfo);
