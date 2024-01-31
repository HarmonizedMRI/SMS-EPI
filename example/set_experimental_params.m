scanner = 'Siemens';  % 'GE' or 'Siemens' (for now!)

if strcmp(scanner, 'GE')
    datfile_suffix = '.h5';    % ScanArchive
else
    datfile_suffix = '.dat';   % Siemens raw data
end

% note trailing slash
%basedir = '~/data/20240124_UM3TUHP_invivo/';
basedir = '~/data/20240130_MGHPrisma_phantom/';
datdir = '';

% SMS-EPI data files
datfile_ghostcal = ['ghostcal' datfile_suffix];    % EPI ghost reference scan
datfile_mb1 = ['mb1' datfile_suffix];              % 2D EPI reference scan for slice GRAPPA
datfile_task = ['task' datfile_suffix];            % 3:30 task fMRI
datfile_rest = ['rest' datfile_suffix];            % 5:13 resting fMRI

% b0 mapping
pfile_b0 = 'P,b0.7';        % 3D GRE dual-echo scan for B0 mapping
deltaTE = 1000/440*1e-3;   % sec

% SMS-EPI acquisition parameters
% The number of temporal frames is not set here; 
% it is determined by the 'runs' parameter on the console
voxelSize = [2.4 2.4 2.4]*1e-3;   % m
nx = 90; ny = 90; nz = 60;        % matrix size
TE = 30e-3;                       % sec
alpha = 52;                       % flip angle (deg)
mb = 6;                           % multiband/SMS factor
pf_ky = 72/90;                    % partial Fourier factor
TR = 0.8;                       % volume TR (sec)

sysGE = toppe.systemspecs();

etl = round(ny*pf_ky);
fov = voxelSize.*[nx ny nz];
np = nz/mb;  % number of partitions (SMS slice groups, or RF shots per temporal frame)
assert(~mod(np,1), 'nz must be a multiple of mb');

% odd/even echo k-space sampling locations (ramp sampling)
[rf,gx,gy,gz,desc,paramsint16,paramsfloat,hdr] = toppe.readmod([basedir 'scanfiles/module7.mod']);
nFID = hdr.rfres;                             % number of data samples in ADC window
del = -1.5;  % apply odd/even k-space delay (samples)
[kxo, kxe] = toppe.utils.getk(sysGE, [basedir 'scanfiles/module7.mod'], nFID, del);

