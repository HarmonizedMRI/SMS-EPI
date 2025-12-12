%% Get MATLAB toolboxes

% Pulseq (+mr namespace)
%system('git clone git@github.com:pulseq/pulseq.git');
%addpath pulseq/matlab


addpath ./lib/
addpath ./lib/
%addpath ~/github/jfnielsen/scanLog/SiemensRead/   % Read Siemens .dat files

addpath /home/jfnielse/Programs/orchestra-sdk-2.1-1.matlab/   % toolbox for loading GE raw files

% hmriutils toolbox
system('git clone --branch develop git@github.com:HarmonizedMRI/utils.git');
addpath utils
 
% +toppe toolbox (misc GE-related scripts)
%system('git clone --branch develop git@github.com:HarmonizedMRI/PulCeq.git');
%addpath PulCeq/matlab
system('git clone --branch develop git@github.com:toppeMRI/toppe.git');
addpath toppe

% 'im' function for display
system('git clone git@github.com:JeffFessler/mirt.git');
cd mirt; 
setup; 
%ir_mex_build_mri;
%ir_mex_build_table;
cd ..;
