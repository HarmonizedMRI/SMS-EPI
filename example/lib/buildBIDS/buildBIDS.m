
srcRoot  = '~/temp';
bidsRoot = '~/dropbox_team/Shared/Data';

sessions = findSessions(srcRoot);

for iSession = 1:numel(sessions)
    info = parseSessionName(sessions(iSession).name, srcRoot, bidsRoot);

    copyBOLD(info);
end
