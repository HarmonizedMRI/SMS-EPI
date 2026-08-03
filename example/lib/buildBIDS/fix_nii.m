function fix_nii(pulseq_filename, product_filename, 
%
% Load a .nii file containing image time-series reconstructed with recon_timeseries.m,
% and write to a new .nii file after updating the following:
%    copy header from product scan, then:
%    set pixeldim 
%    set Datatype = int16
%    scale images to fixed value
%    flip first image dimension    

pulseq_hdr_info = niftiinfo(pulseq_filename);
product_hdr_info = niftiinfo(product_filename);

filename = pulseq_hdr_info.Filename;
filesize = pulseq_hdr_info.Filesize;
%datatype = pulseq_hdr_info.Datatype;
imagesize = pulseq_hdr_info.ImageSize;
%pixeldim = product_hdr_info.PixelDimensions;
pixeldim = [2.4 2.4 2.4 .8];
%transform = pulseq_hdr_info

info = product_hdr_info;
info.Filename = filename;
info.Filesize = filesize;
%info.Datatype = datatype;
info.Datatype = 'int16';
info.ImageSize = imagesize;
info.PixelDimensions = pixeldim;
%info.Transform = transform;
%info.Transform.T(3,3) = 2.4;

data = niftiread(pulseq_hdr_info);
data = data/max(data(:))*2^13;
%data = data/10;
data = int16(data);
data = flipdim(data,1);

outfile = 'test.nii';
niftiwrite(data,outfile,info);
