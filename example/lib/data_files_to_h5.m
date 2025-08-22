% Load ScanArchive/.dat files and write to .h5 files for further processing

% get file names and path
set_experimental_parameters;

% ghost calibration data
D = readraw(datafile_ghostcal, scanner);
hmriutils.epi.io.draw2hdf(D, etl, np, 'ghostcal.h5');

% slice GRAPPA calibration ('ACS') data
D = readraw(datafile_mb1, scanner);
hmriutils.epi.io.draw2hdf(D, etl, np*mb, 'mb1.h5');

% fMRI resting run
D = readraw(datafile_rest, scanner);
% only keep a few frames for testing
%nfr = size(D,3)/etl/np;
%nFramesToKeep = 4;
%hmriutils.epi.io.draw2hdf(D(:,:,1:nFramesToKeep*etl*np), etl, np, ofn); %, 'maxFramesPerFile', 100);
hmriutils.epi.io.draw2hdf(D, etl, np, 'rest.h5', 'maxFramesPerFile', 100);

% fMRI task run(s)
for ii = 1:length(datafile_task)
    D = readraw(datafile_task{ii}, scanner);
    ofn = sprintf('task_%d.h5', ii);
    hmriutils.epi.io.draw2hdf(D, etl, np, ofn, 'maxFramesPerFile', 100);
end
