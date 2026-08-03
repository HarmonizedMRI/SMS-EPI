function ensureDirectory(directoryName)
%ENSUREDIRECTORY Create a directory if it does not already exist.

arguments
    directoryName (1,:) char
end

if isfolder(directoryName)
    return
end

[success, message] = mkdir(directoryName);

if ~success
    error('ensureDirectory:CreateFailed', ...
        'Could not create directory:\n  %s\n%s', ...
        directoryName, ...
        message);
end

end
