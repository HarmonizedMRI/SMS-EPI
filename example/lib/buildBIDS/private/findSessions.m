function sessions = findSessions(srcRoot)
%FINDSESSIONS Return valid source session directories under srcRoot.

arguments
    srcRoot (1,:) char
end

if ~isfolder(srcRoot)
    error('findSessions:MissingSourceRoot', ...
        'Source root does not exist: %s', srcRoot);
end

entries = dir(fullfile(srcRoot, 'sub*'));
entries = entries([entries.isdir]);

names = {entries.name};
keep = ~ismember(names, {'.', '..'});

sessions = entries(keep);

end
