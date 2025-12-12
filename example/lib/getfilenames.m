function F = getfilenames(ifn, vendor)
%
% Get raw data file names from scans.txt file

P = sub_getpaths(ifn);

F = P;

% For GE, scans.txt lists Series folder name only, so get file names here
if strcmp(vendor, 'GE')
    if isfield(P, 'b0')
        D = dir([P.datadir P.b0.name]);
        F.b0.name = [P.b0.name '/' D(end).name];
    end
    if isfield(P, 'epical')
        D = dir([P.datadir P.epical.name]);
        F.epical.name = [P.epical.name '/' D(end).name];
    end
    if isfield(P, 'smscal')
        D = dir([P.datadir P.smscal.name]);
        F.smscal.name = [P.smscal.name '/' D(end).name];
    end
    if isfield(P,'grappacal')
        D = dir([P.datadir P.grappacal.name]);
        F.grappacal.name = [P.grappacal.name '/' D(end).name];
    end
    if isfield(P, 'noise')
        D = dir([P.datadir P.noise.name]);
        F.noise.name = [P.noise.name '/' D(end).name];
    end
    if isfield(P, 'rest')
        for ii = 1:length(P.rest)
            D = dir([P.datadir P.rest(ii).name]);
            F.rest(ii).name = [P.rest(ii).name '/' D(end).name];
        end
    end
    if isfield(P, 'task')
        for ii = 1:length(P.task)
            D = dir([P.datadir P.task(ii).name]);
            F.task(ii).name = [P.task(ii).name '/' D(end).name];
        end
    end
end

return

% read scans.txt file
function P = sub_getpaths(ifn)

L = readlines(ifn);

% first line contains data directory
P.datadir= ['/' strip(L{1}, '/') '/'];  % enforce leading and trailing slash

cnt.rest = 0;
cnt.task = 0;

ii = 2;
while ii <= length(L) & ~isempty(L{ii})
    l = split(L{ii});
    switch(l{1})
        case 'b0'
            P.b0.name = l{2};
        case 'cal'
            P.epical.name = l{2};   % EPI ghost calibration scan (and GE receive gain)
        case 'noise'
            P.noise.name = l{2};
        case '2d'
            P.smscal.name = l{2};   % 2D EPI reference scan for slice GRAPPA recon
        case 'grappacal'
            P.grappacal.name = l{2}; % in-plane grappa reference scan
        case 'task'
            cnt.task = cnt.task + 1;
            P.task(cnt.task).name = l{2};
        case 'rest'
            cnt.rest = cnt.rest + 1;
            P.rest(cnt.rest).name = l{2};
    end
    ii = ii + 1;
end

return
