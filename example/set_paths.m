%% Get MATLAB toolboxes

% Pulseq (+mr namespace)
system('git clone git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab

% hmriutils toolbox
system('git clone git@github.com:HarmonizedMRI/utils.git');
 
% tools for converting to GE scan files
system('git clone git@github.com:HarmonizedMRI/PulCeq.git');
addpath PulCeq/matlab
system('git clone git@github.com:toppeMRI/toppe.git');
addpath toppe
 
