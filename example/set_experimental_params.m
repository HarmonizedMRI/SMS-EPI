
% data files
datdir = '~/mydata/';       % data directory
pfile_cal = 'P,cal.7';      % EPI ghost reference scan (frame 1)
pfile_mb1 = 'P,mb1.7';      % single-slice EPI reference scan for slice GRAPPA
pfile_mb6 = 'P,mb6.7';      % fMRI scan

% SMS-EPI acquisition parameters
voxelSize = [2.4 2.4 2.4]*1e-3;   % m
nx = 90; ny = 90; nz = 60;        % matrix size
TE = 30e-3;                       % sec
alpha = 52;                       % flip angle (deg)
mb = 6;                           % multiband/SMS factor
pf_ky = 1.0;                      % partial Fourier factor
TR = 2*0.8;                       % volume TR (sec)

sysGE = toppe.systemspecs();

etl = round(ny*pf_ky);
fov = voxelSize.*[nx ny nz];
np = nz/mb;  % number of partitions (SMS slice groups, or RF shots per temporal frame)
assert(~mod(np,1), 'nz must be a multiple of mb');

% odd/even echo k-space sampling locations (ramp sampling)
[rf,gx,gy,gz,desc,paramsint16,paramsfloat,hdr] = toppe.readmod([datdir 'module7.mod']);
nFID = hdr.rfres;                             % number of data samples in ADC window
[kxo, kxe] = toppe.utils.getk(sysGE, [datdir 'module7.mod'], nFID);

