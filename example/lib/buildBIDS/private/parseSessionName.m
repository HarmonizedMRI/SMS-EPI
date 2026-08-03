function info = parseSessionName(sessionName, srcRoot, bidsRoot)
%PARSESESSIONNAME Parse an internal session identifier into BIDS paths.
%
% info = parseSessionName(sessionName, srcRoot, bidsRoot)
%
% Example input:
%   sessionName = 'sub00012-umich-750MR-20250115-1'
%
% Example output labels:
%   info.sub = 'sub-00012'
%   info.ses = 'ses-sub00012umich750MR202501151'

arguments
    sessionName (1,:) char
    srcRoot     (1,:) char
    bidsRoot    (1,:) char
end

tokens = regexp( ...
    sessionName, ...
    '^sub(\d+)-(.+)$', ...
    'tokens', ...
    'once');

if isempty(tokens)
    error('parseSessionName:InvalidSessionName', ...
        'Session name does not match the expected format: %s', ...
        sessionName);
end

subjectNumber = str2double(tokens{1});

if isnan(subjectNumber)
    error('parseSessionName:InvalidSubjectNumber', ...
        'Could not parse subject number from: %s', ...
        sessionName);
end

info.sessionName = sessionName;
info.sub = sprintf('sub-%05d', subjectNumber);

sessionValue = regexprep(sessionName, '[^A-Za-z0-9]', '');
info.ses = ['ses-' sessionValue];

info.srcDir = fullfile(srcRoot, sessionName);
info.subjectDir = fullfile(bidsRoot, info.sub);
info.sessionDir = fullfile(info.subjectDir, info.ses);

end
