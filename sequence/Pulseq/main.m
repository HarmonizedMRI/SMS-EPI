% Create SMS EPI scan

% get MATLAB toolboxes
if false
system('git clone git@github.com:pulseq/pulseq.git');
addpath pulseq/matlab
system('git clone git@github.com:HarmonizedMRI/PulCeq.git');
addpath PulCeq/matlab
system('git clone git@github.com:toppeMRI/toppe.git');
addpath toppe
end

% create smsepi.seq
writeSMSEPI;   

% convert to .tar file (smsepi.tar)
smsepi2ge;     

% Plot the .tar file.
% toppe.plotseq() loads the scan files from the MATLA working directory,
% so first we must untar the files.
system('tar xf smsepi.tar');  
toppe.plotseq(sysGE, 'timeRange', [0 0.06]);  % time range in sec
