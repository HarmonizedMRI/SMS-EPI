% Reconstruct Pulseq SMS-EPI fMRI data
% Updated Nov 2024 for new recon workflow

% 1. Add the new exam to sessions.txt
%
% 2. Create the file
%    ../sessions/<subject ID>/<site ID>/<scanner ID>/<date>/<session number>/pulseq/scans.txt
%    listing the raw data file path + names to reconstruct.
%    The folder names must match your new entry in sessions.txt

% get data file names
d = [E.subject '-' E.site '-' E.scanner '-' E.date '-' E.session '/'];
F = getfilenames(['../sessions/' d 'pulseq/scans.txt'], E.vendor);

% set output (.nii file) directory
outputdir = ['./recon-complete/' d 'pulseq/'];
system(sprintf('mkdir -p %s', outputdir));

% odd/even k-space delay (samples) to avoid phase wrap in fit
kspace_delay = E.kspace_delay;

% readout gradient waveform for one EPI echo
readout_trajectory_file = E.readout_trajectory_file;

% SMS-EPI acquisition parameters
% The number of temporal frames is not set here; 
% it is determined by the 'runs' parameter on the console
voxelSize = [3,3,3]*1e-3;   % m
nx = 80; ny = 80; nz = 44;        % matrix size
TE = 30e-3;                       % sec
alpha = 60;                       % flip angle (deg)
mb = 4;                           % multiband/SMS factor
pf_ky = 0.85;                    % partial Fourier factor
TR = 0.9;                       % volume TR (sec)
nTE = 3; % number of echos

%etl = round(ny*pf_ky);
fov = voxelSize.*[nx ny nz];
np = nz/mb;  % number of partitions (SMS slice groups, or RF shots per temporal frame)
assert(~mod(np,1), 'nz must be a multiple of mb');

% odd/even echo k-space sampling locations (ramp sampling)
%%%%%%%%%%%%%%%%%%%%%% tv6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [rf,gx,gy,gz,desc,paramsint16,paramsfloat,hdr] = toppe.readmod(readout_trajectory_file);
% nFID = hdr.rfres;   % number of data samples in ADC window
% sysGE = toppe.systemspecs();
% [kxo, kxe] = toppe.utils.getk(sysGE, readout_trajectory_file, nFID, kspace_delay);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% since pge2, save kxo.kxe into MAT-file using pulseq.calculateKspacePP()
% instead of using readout.mod to get kxo/kxe
load(E.readout_trajectory_file);
nfid = numel(kxo);
kxo = interp1(1:nfid, kxo, (1:nfid) + kspace_delay,'linear','extrap');
kxe = interp1(1:nfid, kxe, (1:nfid) + kspace_delay, 'linear','extrap'); % apply ad-hoc constant delay to avoid phase wrap in the following estimation of a

% b0 mapping parameters(3D GRE)
b0.deltaTE = 1000/440*1e-3;   % sec
b0.N = [100 100 100];         % matrix size
b0.fov = [24 24 24]*1e-2;     % m
b0.nzDummy = 1;               % dummy z loops (for setting receive gain on GE)
%D = dir([datadir '*pulseq_B0.dat']);
%datafile_b0 = [datadir D.name];


