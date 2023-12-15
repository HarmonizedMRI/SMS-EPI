set_paths;

set_experimental_params;

fprintf('Get calibration data...');
get_calibration_data;
fprintf(' done\n');

fprintf('Reconstructing one frame...');
recon_timeseries;
fprintf(' done\n');
