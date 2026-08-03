function copyFileChecked(sourceFilename, destinationFilename)
%COPYFILECHECKED Copy a file and report failures clearly.

arguments
    sourceFilename      (1,:) char
    destinationFilename (1,:) char
end

if ~isfile(sourceFilename)
    error('copyFileChecked:MissingSource', ...
        'Source file does not exist: %s', sourceFilename);
end

ensureDirectory(fileparts(destinationFilename));

[success, message] = copyfile( ...
    sourceFilename, ...
    destinationFilename, ...
    'f');

if ~success
    error('copyFileChecked:CopyFailed', ...
        'Could not copy:\n  %s\nTo:\n  %s\n%s', ...
        sourceFilename, ...
        destinationFilename, ...
        message);
end

fprintf('Copied:\n  %s\n', destinationFilename);

end
