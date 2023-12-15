
% data files
datdir = '/mnt/storage/jfnielse/HarmonizedMRI/EPI-test/20231204_UM3TMR750_phantom/';
pfileSensGRE = 'P,b0.7';    % 3D GRE for B0 field mapping and (potentially) ACS data for recon 
pfileCal = 'P,cal.7';       % first frame is EPI reference scan
pfile_mb6   = 'P,mb6.7';
pfile_mb1   = 'P,2d.7';     % single-slice EPI reference scan for slice GRAPPA

% acquisition parameters
res = [2.4 2.4 2.4]*1e-3;   % voxel dimensions (m)
N = [90 90 60];             % image matrix size
etl = 90;                   % no partial Fourier (for now)
mb = 6;                     % multiband/SMS factor

% Same readout for all smsepi scans
readoutFile = [datdir 'scanfiles/module7.mod'];
%caipiFile = [datdir 'scanfiles/caipi.mat'];                  % created by skippedcaipi_sampling.py

% tool paths
addpath ~/github/HarmonizedMRI/utils/ 

sysGE = toppe.systemspecs();

fov = res.*N;
fovXcm = fov(1)*100;
[nx ny nz] = deal(N(1), N(2), N(3));

np = nz/mb;  % number of partitions (SMS slice groups, or RF shots per frame)

