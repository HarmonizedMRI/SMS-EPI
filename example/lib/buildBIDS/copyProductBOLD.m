function copyProductBOLD(srcDir, bidsDir)

% Convert product NIfTI files into BIDS layout.
%
% Example:
%
% copyProductBOLD( ...
%     '/data/internal', ...
%     '/data/BIDS');
%
% Expected source layout:
%
%   srcDir/
%       sub00012-umich-750MR-20250115-1/
%           product/
%               task_run1.h5.nii
%               task_run2.h5.nii
%
% Output:
%
%   bidsDir/
%       sub-00012/
%           ses-sub00012umich750MR202501151/
%               func/
%                   sub-00012_ses-sub00012umich750MR202501151_task-rest_acq-product_run-01_bold.nii

%% find all session folders

D = dir(fullfile(srcDir,'sub*'));

for i = 1:numel(D)

    if ~D(i).isdir
        continue
    end

    sessionName = D(i).name;

    %% parse folder name

    tok = regexp(sessionName,...
        '^sub(\d+)-(.+)$',...
        'tokens','once');

    if isempty(tok)
        warning('Skipping %s',sessionName);
        continue
    end

    subjectNumber = tok{1};

    subLabel = sprintf('sub-%04d',str2double(subjectNumber));

    % Remove hyphens to make valid BIDS label
    sesLabel = ['ses-' regexprep(sessionName,'-','')];

    %% source folder

    srcFunc = fullfile(srcDir,sessionName,'product');

    if ~isfolder(srcFunc)
        continue
    end

    %% destination

    dstFunc = fullfile( ...
        bidsDir,...
        subLabel,...
        sesLabel,...
        'func');

    if ~exist(dstFunc,'dir')
        mkdir(dstFunc);
    end

    %% copy all task files

    F = dir(fullfile(srcFunc,'task_run*.h5.nii'));

    for j = 1:numel(F)

        runTok = regexp(F(j).name,...
            'run(\d+)',...
            'tokens','once');

        run = str2double(runTok{1});

        outName = sprintf('%s_%s_task-rest_acq-product_run-%02d_bold.nii',...
            subLabel,...
            sesLabel,...
            run);

        copyfile( ...
            fullfile(srcFunc,F(j).name), ...
            fullfile(dstFunc,outName));

        fprintf('Copied %s\n',outName);

    end

end
