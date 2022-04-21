function epiInfo = caipiepifmri(scanType, N, FOV, flip, nFrames, varargin)
% function epiInfo = caipiepifmri(scanType, N, FOV, flip, nFrames, varargin)
%
% Creat 3D/SMS CAIPI EPI fMRI scan in TOPPE
%
% 3D version implements sequence in Narsude et al MRM 2016, 
% "Three-Dimensional Echo Planar Imaging with Controlled Aliasing:
% A Sequence for High Temporal Resolution Functional MRI"
%
% SMS version is similar except partition-encode replaced
% by z shift of excited slices.
%
% Inputs:
%  scanType   string   '3D' or 'SMS'
%  N          [1 3]    matrix size
%  FOV        [1 3]    field of view (cm)
%  flip       [1 1]    flip angle (degrees)
%  nFrames    [1 1]    number of time frames
%
% Optional keywords arguments with defaults:
% 
%
% Outputs:
%  epiInfo    struct with various information needed for reconstruction
% 
% In addition, this script creates the file 'fmri.tar', that contains
% the following scan files:
%      toppeN.entry
%      seqstamp.txt
%      scanloop.txt
%      modules.txt
%      .mod files
% These files are described here: https://github.com/toppeMRI/toppe/blob/main/Files.md
% These files are also written to the current Matlab working folder,
% which allows you to plot the scan right away.
%
% Example usage:
%   >>  = caipiepifmri('SMS');
%   >> toppe.playseq(4,sys, 'tpause', 0.2);

%% Set defaults for keyword arguments and replace if provided
arg.entryFile = 'toppeN.entry';
arg.filePath = '/usr/g/research/pulseq/fmri/';  % location of scan files on scanner host
arg.rfSpoilSeed = 117;         % RF spoiling phase increment factor (degrees)
arg.mb = 4;        % sms/multiband factor (number of simultaneous slices)
arg.Ry = 1;        % y acceleration factor
arg.pf_ky = 1.0;   % Partial Fourier factor
arg.Rz = arg.mb;      % z acceleration factor
arg.Delta = arg.mb;   % kz step size (multiples of 1/fov(3))
arg.nCalFrames = 4;   % Odd/even echo calibration frames (ky and kz-encoding off)
arg.doSim = true;     % Plot simulated SMS slice profile
arg.alternateReadout = false;  % flip sign of gx gradient every time-frame
arg.nSpoilCycles = 2;   % number of cycles/voxel of gradient spoiling along z

% Hardware limits (for design and detailed timing calculations).
% 'maxSlew' and 'maxGrad' options can be < scanner limit, 
% and can vary across .mod files. 
arg.sys = toppe.systemspecs('maxSlew', 20, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
    'myrfdel', 200, ...  % psd_rf_wait
    'daqdel', 100, ...   % psd_grd_wait
    'timessi', 100, ...
    'gradient', 'hrmb');  

% substitute with provided keyword arguments
arg = toppe.utils.vararg_pair(arg, varargin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EDIT THIS SECTION AS NEEDED

% TODO: specify TR and TE

% slab/slice excitation pulse parameters
if N(3) > 1  % 3D
    seq.slab.thick = 0.8*FOV(3); % to avoid wrap-around in z
    seq.slab.tbw = 8;
    seq.slab.type = 'st';  % 'st' = small-tip. 'ex' = 90 degree design
    seq.slab.ftype = 'min'; % minimum-phase SLR pulse is well suited for 3D slab excitation
    seq.slab.dur = 1.5;  % ms
else
    seq.slab.thick = FOV(3); 
    seq.slab.tbw = 8;
    seq.slab.type = 'st';  % 'st' = small-tip. 'ex' = 90 degree design
    seq.slab.ftype = 'ls'; % minimum-phase SLR pulse is well suited for 3D slab excitation
    seq.slab.dur = 1.5;  % ms
end

% SMS excitation pulse parameters
seq.sms.type = 'st';      % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
seq.sms.ftype = 'ls';     
seq.sms.tbw = 6;          % time-bandwidth product
seq.sms.dur = 4;          % msec
seq.sms.sliceSep = FOV(3)/arg.mb;   % center-to-center separation between SMS slices (cm)

% slew rates for waveform design (G/cm/ms)
slewRead = [11 15 15];  % x/y/z max slew for EPI train
slewPre = 8;            % for prewinder trapezoid

%% DONE WITH CUSTOM EDITS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmp(scanType, '3D') & ~strcmp(scanType, 'SMS')
    error('Supported scan types: 3D, SMS');
end

seq.voxelSize = FOV./N;

nx = N(1);
ny = N(2);
nz = N(3);

if mod(ny,arg.Ry) > 0
    error('ny/Ry must be an integer');
end
if mod(nz,arg.Rz) > 0
    error('nz/Rz must be an integer');
end


%% Create modules.txt and entryFile

% Entries are tab-separated
fid = fopen('modules.txt', 'wt');
fprintf(fid, 'Total number of unique cores\n');
fprintf(fid, '%d\n', 3);
fprintf(fid, 'fname  duration(us)    hasRF?  hasDAQ?\n');
fprintf(fid, '%s\t0\t1\t0\n', 'tipdown.mod');
fprintf(fid, '%s\t0\t0\t0\n', 'prephase.mod');
fprintf(fid, '%s\t0\t0\t1\n', 'readout.mod');
fclose(fid);

% Write entry file.
% This can be edited by hand as needed after copying to scanner.
toppe.writeentryfile(arg.entryFile, ...
    'filePath', arg.filePath);


%% Slab/SMS excitation
if strcmp(scanType, 'SMS')
    tmp = arg.sys;
    tmp.maxSlew = 8;   % G/cm/ms. Reduce PNS during slice rephaser.
    [ex.rf, ex.g, freq] = getsmspulse(flip, seq.voxelSize(3), seq.sms.tbw, seq.sms.dur, ...
        arg.mb, seq.sms.sliceSep, tmp, ...
	    'doSim', arg.doSim, ...   % Plot simulated SMS slice profile
	    'type', seq.sms.type, ...
	    'ftype', seq.sms.ftype);

    toppe.writemod(arg.sys, 'rf', ex.rf, 'gz', ex.g, ...
	    'ofname', 'tipdown.mod' );

    freq = freq/seq.sms.sliceSep*seq.voxelSize(3); % frequency offset for z shift of seq.voxelSize(3)
else
    % slab select
    toppe.utils.rf.makeslr(flip, seq.slab.thick, ...
        seq.slab.tbw, seq.slab.dur, nz*arg.nSpoilCycles, arg.sys, ...
        'type', seq.slab.type, ...     % 'st' = small-tip. 'ex' = 90 degree design
        'ftype', seq.slab.ftype, ...  
        'spoilDerate', 0.5, ...
        'ofname', 'tipdown.mod');
end


%% Readout 
[gx, gy, gz, gpre, esp, gx1, kz, nBlipMax] = getcaipiepireadout(FOV, N, ...
    arg.Ry, arg.pf_ky, ...
    arg.Rz, arg.Delta, ...
    arg.sys.maxGrad, slewRead, slewPre, ...
    arg.sys.raster*1e3, arg.sys.forbiddenEspRange);

toppe.writemod(arg.sys, ...
    'gx', gpre.x, 'gy', gpre.y, 'gz', gpre.z, ...
    'ofname', 'prephase.mod' );

toppe.writemod(arg.sys, ...
    'gx', gx, 'gy', gy, 'gz', gz, ...
    'ofname', 'readout.mod' );


%% Scan loop
rfphs = 0;              % radians
rfSpoilSeed_cnt = 0;

toppe.write2loop('setup', arg.sys, 'version', 4);   % Initialize scanloop.txt

IZ = 1:arg.Rz:nz;

% temporal loop
for ifr = 1:(arg.nCalFrames  + nFrames)
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bFrame %d of %d', ifr, arg.nCalFrames + nFrames);
    isCalScan = (ifr < arg.nCalFrames + 1);

    % shift ky sampling pattern (for Ry = 2, this corresponds to UNFOLD)
    %iy = (1 + (-1)^(ifr))/2 + 1;   % 1 or 2
    %a_gy = ((iy-1+0.5)-ny/2)/(ny/2);
    a_gy = (1-isCalScan)*((1-1+0.5)-ny/2)/(ny/2);

    if ~isCalScan & arg.alternateReadout
        a_gx = (-1)^(ifr+1); % alternate so we can combine frames to simulate flyback EPI
    else
        a_gx = 1;
    end

    % z encoding / SMS slice shift loop
    for ii = 1:length(IZ)  

        if strcmp(scanType, 'SMS')
		    f = round((ii-0.5-length(IZ)/2)*freq);  % frequency offset (Hz) for slice shift
            a_gz = 0;
        else
            f = 0;
            a_gz = (1-isCalScan)*((IZ(ii)-1+0.5)-nz/2)/(nz/2);
        end

        % rf excitation
        toppe.write2loop('tipdown.mod', arg.sys, ...
            'RFoffset', f, ...
            'RFphase', rfphs);

        % readout 
        % data is stored in 'slice', 'echo', and 'view' indeces
        toppe.write2loop('prephase.mod', arg.sys, ...
            'Gamplitude', [-1.0*a_gx a_gy a_gz]');
        toppe.write2loop('readout.mod', arg.sys, ...
            'Gamplitude', [a_gx (1-isCalScan) (1-isCalScan)]', ...
            'DAQphase', rfphs, ...
            'slice', ii, 'view', ifr);

        % make net gradient area per TR constant (to achieve steady state)
        toppe.write2loop('prephase.mod', arg.sys, ...
            'Gamplitude', [0 -a_gy -a_gz]');

        % update rf phase (RF spoiling)
        rfphs = rfphs + (arg.rfSpoilSeed/180 * pi)*rfSpoilSeed_cnt ;  % radians
        rfSpoilSeed_cnt = rfSpoilSeed_cnt + 1;
    end
end
toppe.write2loop('finish', arg.sys);
fprintf('\n');

%% Create 'sequence stamp' file for TOPPE.
% This file is listed in the 5th row in entryFile
% NB! The file entryFile must exist in the folder from where this script is called.
toppe.preflightcheck(arg.entryFile, 'seqstamp.txt', arg.sys);

%% create tar file
system(sprintf('tar cf fmri.tar %s seqstamp.txt scanloop.txt modules.txt *.mod gx1.mat kz.mat', arg.entryFile));

toppe.utils.scanmsg(arg.entryFile);


%% Return struct needed for recon
epiInfo.N = N;
epiInfo.FOV = FOV;
epiInfo.gpre = gpre;  % x prewinder (G/cm)
epiInfo.gx = gx;      % full echo train (G/cm)
epiInfo.gx1 = gx1;    % one echo
epiInfo.nBlipMax = nBlipMax;  % number of samples durings turns

