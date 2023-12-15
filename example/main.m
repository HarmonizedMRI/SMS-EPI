% get the code
set_paths;

% Edit the following script as needed
set_experimental_params;

% Create .seq files
write_sequences;

% After scanning, place the raw data files in the folder 
% specified in set_experimental_params.m.
% Then run the remainder of this script.

% Get calibration data and reconstruct 
get_calibration_data;
recon_timeseries;
