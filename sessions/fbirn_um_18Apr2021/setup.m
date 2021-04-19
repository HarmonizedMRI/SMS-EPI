addpath ~/gitlab/fMRI/toppe/
curdir = pwd; cd ~/github/mirt; setup; cd(curdir);

% Add the following to ~/.bashrc:
% TOOLBOX_PATH=/home/jon/Programs/bart-0.7.00
addpath(strcat(getenv('TOOLBOX_PATH'), '/matlab'));
