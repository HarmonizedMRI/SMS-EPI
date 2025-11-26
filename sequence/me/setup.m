% Pulseq +mr toolbox
system('git clone --branch v1.5.1 git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab

% Code needed to convert to GE sequence representation
%system('git clone --branch tv7 git@github.com:HarmonizedMRI/PulCeq.git');
%addpath PulCeq/matlab
%addpath PulCeq/matlab/DataHash
addpath ~/github/HarmonizedMRI/PulCeq/matlab
addpath ~/github/HarmonizedMRI/PulCeq/matlab/DataHash

% toppe.utils.makeslr
system('git clone --branch v1.9.1 git@github.com:toppeMRI/toppe.git');
addpath toppe

% Python code for creating CAIPI ky-kz sampling pattern
system('git clone git@github.com:HarmonizedMRI/3DEPI.git');

