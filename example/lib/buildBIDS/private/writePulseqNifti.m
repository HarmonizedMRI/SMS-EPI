function writePulseqNifti(pulseqFilename, productFilename, outputFilename, overwrite)
%WRITEPULSEQNIFTI Apply product geometry to a Pulseq reconstruction.
%
% The Pulseq data are:
%   - flipped along dimension 1;
%   - globally scaled to maximum absolute magnitude 2^13;
%   - converted to int16.
%   - pixelDimensions are set to [2.4 2.4 2.4 0.8]
%
% The output header is based on the corresponding product NIfTI header.

arguments
    pulseqFilename  (1,:) char
    productFilename (1,:) char
    outputFilename  (1,:) char
    overwrite       (1,1) logical = false
end

if nargin < 3
    overwrite = false;
end

if isfile(outputFilename)
    if overwrite
        fprintf('Overwriting:\n  %s\n', outputFilename);
    else
        fprintf('Skipping existing file:\n  %s\n', outputFilename);
        return
    end
end


if ~isfile(pulseqFilename)
    error('writePulseqNifti:MissingPulseqFile', ...
        'Pulseq file does not exist: %s', pulseqFilename);
end

if ~isfile(productFilename)
    error('writePulseqNifti:MissingProductFile', ...
        'Product reference file does not exist: %s', productFilename);
end

pulseqInfo = niftiinfo(pulseqFilename);
productInfo = niftiinfo(productFilename);

data = niftiread(pulseqInfo);

validateSpatialDimensions( ...
    size(data), ...
    productInfo.ImageSize, ...
    pulseqFilename, ...
    productFilename);

data = flip(data, 1);
data = scaleToInt16(data, 2^13);

outputInfo = productInfo;
outputInfo.Datatype = 'int16';
outputInfo.ImageSize = size(data);

if ndims(data) >= 4
    outputInfo.PixelDimensions = [2.4 2.4 2.4 0.8];
else
    outputInfo.PixelDimensions = [2.4 2.4 2.4];
end

% These fields describe the source file and should not be reused.
if isfield(outputInfo, 'Filename')
    outputInfo.Filename = '';
end

if isfield(outputInfo, 'Filesize')
    outputInfo.Filesize = 0;
end

ensureDirectory(fileparts(outputFilename));

niftiwrite( ...
    data, ...
    outputFilename, ...
    outputInfo, ...
    'Compressed', false);

fprintf('Created corrected Pulseq NIfTI:\n  %s\n', outputFilename);

end
