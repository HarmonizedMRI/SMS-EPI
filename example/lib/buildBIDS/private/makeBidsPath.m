function [filename, destination] = makeBidsPath(info, modality, suffix)
%MAKEBIDSPATH Construct a BIDS filename and destination path.
%
% [filename, destination] = makeBidsPath(info, modality, suffix)
%
% Example:
%   suffix = 'task-rest_acq-product_run-01_bold.nii'
%
% Output:
%   sub-00012_<session>_task-rest_acq-product_run-01_bold.nii

arguments
    info     (1,1) struct
    modality (1,:) char
    suffix   (1,:) char
end

validateSessionInfo(info);

if isempty(modality)
    error('makeBidsPath:EmptyModality', ...
        'Modality directory cannot be empty.');
end

if isempty(suffix)
    error('makeBidsPath:EmptySuffix', ...
        'BIDS filename suffix cannot be empty.');
end

filename = sprintf('%s_%s_%s', ...
    info.sub, ...
    info.ses, ...
    suffix);

destinationDir = fullfile(info.sessionDir, modality);
ensureDirectory(destinationDir);

destination = fullfile(destinationDir, filename);

end
