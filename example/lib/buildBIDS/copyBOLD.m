function copyBOLD(info)
%COPYBOLD Copy product and Pulseq BOLD images for one session.
%
% copyBOLD(info)
%
% Required fields in info:
%   info.srcDir      Source session directory
%   info.sub         BIDS subject label, e.g. 'sub-00012'
%   info.ses         BIDS session label
%   info.sessionDir  BIDS destination session directory
%
% Expected source layout:
%
%   info.srcDir/
%       product/
%           task_run1.h5.nii
%           task_run2.h5.nii
%       pulseq/
%           task_run1.h5.nii
%           task_run2.h5.nii
%
% Output:
%
%   info.sessionDir/func/
%       sub-00012_<session>_task-rest_acq-product_run-01_bold.nii
%       sub-00012_<session>_task-rest_acq-pulseq_run-01_bold.nii
%
% Product images are copied without modification.
%
% Pulseq images are:
%   - assigned the spatial header from the matching product image;
%   - flipped along image dimension 1;
%   - scaled globally to a maximum absolute value of 2^13;
%   - written as int16;
%   - assigned voxel dimensions [2.4 2.4 2.4] mm and TR 0.8 s.

arguments
    info (1,1) struct
end

validateSessionInfo(info);

productDir = fullfile(info.srcDir, 'product');
pulseqDir  = fullfile(info.srcDir, 'pulseq');

if ~isfolder(productDir)
    warning('copyBOLD:MissingProductDir', ...
        'Product directory not found: %s', productDir);
    return
end

funcDir = fullfile(info.sessionDir, 'func');
ensureDirectory(funcDir);

productFiles = dir(fullfile(productDir, 'task_run*.h5.nii'));

if isempty(productFiles)
    warning('copyBOLD:NoProductFiles', ...
        'No product BOLD files found in: %s', productDir);
    return
end

for iFile = 1:numel(productFiles)

    sourceName = productFiles(iFile).name;
    runNumber = parseRunNumber(sourceName);

    if isempty(runNumber)
        warning('copyBOLD:InvalidFilename', ...
            'Skipping unrecognized BOLD filename: %s', sourceName);
        continue
    end

    productSource = fullfile(productDir, sourceName);

    productSuffix = sprintf( ...
        'task-rest_acq-product_run-%02d_bold.nii', ...
        runNumber);

    [~, productDestination] = makeBidsPath( ...
        info, ...
        'func', ...
        productSuffix);

    copyFileChecked(productSource, productDestination);

    pulseqSource = fullfile(pulseqDir, sourceName);

    if ~isfile(pulseqSource)
        warning('copyBOLD:MissingPulseqFile', ...
            'Matching Pulseq file not found: %s', pulseqSource);
        continue
    end

    pulseqSuffix = sprintf( ...
        'task-rest_acq-pulseq_run-%02d_bold.nii', ...
        runNumber);

    [~, pulseqDestination] = makeBidsPath( ...
        info, ...
        'func', ...
        pulseqSuffix);

    writePulseqNifti( ...
        pulseqSource, ...
        productSource, ...
        pulseqDestination);

end

end
