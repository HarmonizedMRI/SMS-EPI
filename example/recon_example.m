
% Get the code
get_code_and_set_paths;

% Get experimental parameters
set_experimental_params;

% Load ScanArchive files and write to custom .h5 file
data_files_to_h5;

% Get EPI ghost calibration data 
get_ghost_calibration_data;

% Get individual slice (ACS) data for slice GRAPPA
get_acs_data;

% Reconstruct fMRI data
recon_timeseries;
