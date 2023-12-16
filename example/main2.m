% Get experimental parameters
set_experimental_params;

% Get EPI ghost calibration data 
get_calibration_data;

% Get individual slice (ACS) data for slice GRAPPA
get_acs_data;

% Reconstruct fMRI data
recon_timeseries;
