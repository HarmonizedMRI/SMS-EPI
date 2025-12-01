% main SMS-EPI fMRI image reconstruction script
%
% 1. Add exam to sessions.txt
% 2. Create the file sessions/subject-site-scanner-date-session/pulseq/scans.txt file listing Pulseq file names
% 3. Run this script:
%    Set booleans at start of script to select actions, and row number in sessions.txt file. 
%    >> rownum = 1;  % row number in sessions.txt file listing scan sessions
%    >> recon_example;
%
% Reconstructed images are saved in recon-complete/subject-site-scanner-date-session/

rownum = 10;  % row number in sessions.txt (list of scan sessions)

auto = false;  % pause after each step before continuing

doGhostCal = 1;
getACS = 1; 

newCal = 1;   % 0: before Nov 2024. 1: combines EPI cal and 2D GRAPPA ref/cal scan

reconTask = 0;
reconRest = 1;
nRestToRecon = 1;

% set paths
%if ~exist('pathsHaveBeenSet'
get_code_and_set_paths; 

% Set exam to reconstruct
E = getexaminfo('sessions.txt', rownum)

% where to put intermediate files
tmpdir = ['./tmp' num2str(rownum) '/'];
system(sprintf('mkdir -p %s', tmpdir));

set_experimental_parameters;

% B0 map
%get_b0;

% EPI ghost calibration. Saves linear ghost correction parameters in a.mat
% If you observe phase wraps in plots, change kspace_delay in getexaminfo.m
if doGhostCal
datafile_ghostcal = [F.datadir F.epical.name];
D = readraw(datafile_ghostcal, E.vendor);
h5file_ghostcal = [tmpdir 'ghostcal.h5'];
hmriutils.epi.io.draw2hdf(D, etl, np, h5file_ghostcal);
get_ghost_calibration_data;  
if ~auto
    input('Hit Enter to continue ');
end
end

% Get slice GRAPPA calibration data (dcal)
if getACS
datafile_mb1 = [F.datadir F.smscal.name];
D = readraw(datafile_mb1, E.vendor);
h5file_mb1 = [tmpdir 'mb1.h5'];
hmriutils.epi.io.draw2hdf(D, etl, np*mb, h5file_mb1);
if newCal %& strcmp(E.vendor, 'Siemens')
get_acs_data_2;
else
get_acs_data;
end
if ~auto
    input('Hit Enter to continue ');
end
end

% Reconstruct fMRI task runs
loadData = true;
if reconTask & isfield(F, 'task')
    for ii = 1 : length(F.task)
        ifn = ['task_run' num2str(ii) '.h5'];   % do not include tmpdir path here
        if loadData
            D = readraw([F.datadir F.task(ii).name], E.vendor);
            hmriutils.epi.io.draw2hdf(D, etl, np, [tmpdir ifn], 'maxFramesPerFile', 50);
        end
        nFrames = 264;
        recon_timeseries;
    end
end

% Reconstruct fMRI rest runs
loadData = true;
if reconRest & isfield(F, 'rest')
    for ii = 1 : min(nRestToRecon, length(F.rest))
        ifn = ['rest_run' num2str(ii) '.h5'];   % do not include tmpdir path here
        if loadData
            D = readraw([F.datadir F.rest(ii).name], E.vendor);
            hmriutils.epi.io.draw2hdf(D, etl, np, [tmpdir ifn], 'maxFramesPerFile', 50);
        end
        nFrames = 2;
        recon_timeseries;
    end
end

return

