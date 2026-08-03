function copyFileChecked(sourceFilename, destinationFilename, overwrite)

if nargin < 3
    overwrite = false;
end

if isfile(destinationFilename)
    if overwrite
        fprintf('Overwriting:\n  %s\n', destinationFilename);
    else
        fprintf('Skipping existing file:\n  %s\n', destinationFilename);
        return
    end
end

ensureDirectory(fileparts(destinationFilename));

[success,msg] = copyfile(sourceFilename,destinationFilename,'f');

if ~success
    error('Could not copy:\n%s',msg);
end

fprintf('Copied:\n  %s\n', destinationFilename);
