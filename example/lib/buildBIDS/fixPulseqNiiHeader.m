function fixPulseqNiiHeader(pulseq_filename, product_filename, outfile)
%FIXPULSEQNIIHEADER Correct the header and orientation of a reconstructed Pulseq NIfTI.
%
% fix_nii(pulseq_filename, product_filename, outfile)
%
% The function:
%   1. Reads image data from pulseq_filename.
%   2. Uses product_filename as the reference for spatial geometry.
%   3. Sets the voxel dimensions to [2.4 2.4 2.4] mm and TR to 0.8 s.
%   4. Flips the first image dimension.
%   5. Scales the data by a fixed factor (10)
%   6. Writes the result as int16.
%
% Example:
%   fix_nii( ...
%       'pulseq/task_run1.nii', ...
%       'product/task_run1.h5.nii', ...
%       'task_run1_fixed.nii');

arguments
    pulseq_filename  (1,:) char
    product_filename (1,:) char
    outfile          (1,:) char
end

%% Read headers and Pulseq image data

pulseq_info = niftiinfo(pulseq_filename);
product_info = niftiinfo(product_filename);

data = niftiread(pulseq_info);

%% Check spatial dimensions

pulseq_size = size(data);
product_size = product_info.ImageSize;

% MATLAB may omit trailing singleton dimensions.
if numel(pulseq_size) < 3
    pulseq_size(end+1:3) = 1;
end

if numel(product_size) < 3
    product_size(end+1:3) = 1;
end

assert(isequal(pulseq_size(1:3), product_size(1:3)), ...
    ['Spatial image dimensions do not match.\n' ...
     'Pulseq:  [%s]\n' ...
     'Product: [%s]'], ...
    num2str(pulseq_size(1:3)), ...
    num2str(product_size(1:3)));

%% Flip reconstructed data into the product-image orientation

data = flip(data, 1);

%% Scale and convert to int16

data = double(data);

scale_target = 2^13;
max_value = max(abs(data(:)));

if max_value > 0
    data = data .* (scale_target / max_value);
else
    warning('Pulseq image contains only zeros; no scaling was applied.');
end

% Ensure values remain in the int16 range.
%data = max(min(data, double(intmax('int16'))), ...
%               double(intmin('int16')));
data = data/10;

data = int16(round(data));

%% Construct output header

% Start with the product header so that orientation and spatial position
% are inherited from the product scan.
info = product_info;

info.Datatype = 'int16';
info.ImageSize = size(data);

% niftiwrite expects one PixelDimensions value per image dimension.
if ndims(data) >= 4
    info.PixelDimensions = [2.4 2.4 2.4 0.8];
else
    info.PixelDimensions = [2.4 2.4 2.4];
end

%% Write corrected NIfTI

% The data have already been flipped to match the product geometry.
% Therefore, retain the product Transform without modifying it.
niftiwrite(data, outfile, info, 'Compressed', false);

fprintf('Wrote corrected NIfTI:\n  %s\n', outfile);

return



function old_fix_nii(pulseq_filename, product_filename) 
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
