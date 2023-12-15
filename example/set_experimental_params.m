
% data files
datdir = '/mnt/storage/jfnielse/HarmonizedMRI/EPI-test/20231204_UM3TMR750_phantom/';
pfileSensGRE = 'P,b0.7';    % 3D GRE for B0 field mapping and (potentially) ACS data for recon 
pfileCal = 'P,cal.7';       % first frame is EPI reference scan
pfile_mb6   = 'P,mb6.7';
pfile_mb1   = 'P,2d.7';     % single-slice EPI reference scan for slice GRAPPA

% SMS-EPI acquisition parameters
voxelSize = [2.4 2.4 2.4]*1e-3;   % m
TE = 30e-3;                       % sec
alpha = 52;                       % flip angle (deg)
nx = 90; ny = 90; nz = 60;        % matrix size
mb = 6;                           % multiband/SMS factor
pf_ky = 1.0;                      % partial Fourier factor
TR = 2*0.8;                       % volume TR (sec)

fov = res.*[nx ny nz];
np = nz/mb;  % number of partitions (SMS slice groups, or RF shots per temporal frame)
assert(~mod(np,1), 'nz must be a multiple of mb');

% odd/even echo k-space sampling locations (ramp sampling)
readoutFile = [datdir 'scanfiles/module7.mod'];
%caipiFile = [datdir 'scanfiles/caipi.mat'];                  % created by skippedcaipi_sampling.py

sysGE = toppe.systemspecs();

fovXcm = fov(1)*100;


