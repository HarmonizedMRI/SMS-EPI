
curdir = pwd; cd ~/github/mirt; setup; cd(curdir)        % MIRT

%% BART setup (to calculate sensitivity maps)
% Add the following to ~/.bashrc:
% TOOLBOX_PATH=/home/jon/Programs/bart-0.7.00
%>> setenv('TOOLBOX_PATH', '/home/jon/Programs/bart-0.7.00'); % if getenv doesn't find it
%>> addpath(strcat(getenv('TOOLBOX_PATH'), '/matlab'));

