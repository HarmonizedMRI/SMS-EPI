% Pulseq +mr toolbox
system('git clone --branch v1.5.1 git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab
warning('OFF', 'mr:restoreShape');

% Code needed to convert to GE sequence representation
system('git clone --branch main git@github.com:HarmonizedMRI/PulCeq.git');
addpath PulCeq/matlab
addpath PulCeq/matlab/DataHash

system('git clone --branch main git@github.com:HarmonizedMRI/utils.git');
addpath utils

% toppe.utils.makeslr
system('git clone --branch v1.9.1 git@github.com:toppeMRI/toppe.git');
addpath toppe

