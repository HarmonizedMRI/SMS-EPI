% Pulseq +mr toolbox
system('git clone --branch v1.5.1 git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab
warning('OFF', 'mr:restoreShape');

% Pulseq on GE toolbox
%system('git clone --branch main git@github.com:HarmonizedMRI/PulCeq.git');
%addpath PulCeq/matlab
%addpath PulCeq/matlab/DataHash
addpath ~/github/HarmonizedMRI/PulCeq/matlab
addpath ~/github/HarmonizedMRI/PulCeq/matlab/DataHash

% function for creating SLR and SMS excitation pulses
system('git clone --branch main git@github.com:HarmonizedMRI/utils.git');
addpath utils
system('git clone --branch v1.9.1 git@github.com:toppeMRI/toppe.git');
addpath toppe

