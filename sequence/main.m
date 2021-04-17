addpath ~/github/pulseq/matlab/
addpath ~/gitlab/fMRI/toppe/
addpath ~/github/toppeMRI/PulseGEq/

% make 3D B0 mapping/coild sensitivity mapping sequence
gre3d;

% make SMS EPI fMRI sequence
fmri2depi;
