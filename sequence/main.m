% paths for JFN
addpath ~/github/pulseq/matlab/
addpath ~/gitlab/fMRI/toppe/
addpath ~/github/toppeMRI/PulseGEq/

% paths for MZ
%addpath ~/pulseq_home/github/pulseq/matlab/
%addpath ~/pulseq_home/github/toppe/
%addpath ~/pulseq_home/github/PulseGEq/

% make 3D B0 mapping/coild sensitivity mapping sequence
gre3d;

% 2D EPI fMRI sequence
fmri2depi;

% make SMS EPI fMRI sequence
smsepi;
