function [seq, sys] = caipiepifmri(scanType)
% 3D/SMS CAIPI EPI fMRI scan in TOPPE
%
% 3D version implements sequence in Narsude et al MRM 2016, 
% "Three-Dimensional Echo Planar Imaging with Controlled Aliasing:
% A Sequence for High Temporal Resolution Functional MRI"
%
% SMS version is similar except partition-encode replaced
% by z shift of excited slices.
%
% Inputs:
%  scanType    '3D' or 'SMS'
%
% Outputs:
%  TOPPE scan files        modules.txt, scanloop.txt, tipdown.mod, and readout.mod

if ~strcmp(scanType, '3D') & ~strcmp(scanType, 'SMS')
    error('Supported scan types: 3D, SMS');
end

%% Set hardware limits (for design and detailed timing calculations).
% 'maxSlew' and 'maxGrad' options can be < scanner limit, and can vary across .mod files. 
sys.ge = toppe.systemspecs('maxSlew', 20, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
    'myrfdel', 200, ...  % psd_rf_wait
    'daqdel', 100, ...   % psd_grd_wait
    'timessi', 100, ...
    'gradient', 'xrm');  % MR750 scanner. Used to plot PNS profile in toppe.plotmod()


%% Sequence parameters

seq.imSize = [80 80 40];
seq.fov = 0.24*seq.imSize;    % cm
seq.voxelSize = seq.fov./seq.imSize;

seq.nFrames = 20;  % number of fMRI frames
seq.nCalFrames = 4;  % Odd/even echo calibration frames (ky and kz-encoding turned off)

seq.nSpoilCycles = 2;   % number of cycles of gradient spoiling along z

% TODO: specify TR and TE

% EPI CAIPI sampling parameters
mb = 4;      % sms/multiband factor (number of simultaneous slices)
seq.Ry = 1;       % acceleration factor
seq.pf_ky = 1.0;  % Partial Fourier factor
seq.Rz = mb;
seq.Delta = mb;  % kz step size (multiples of 1/fov(3))

% slab excitation pulse parameters
seq.slab.flip = 15;  % flip angle (degrees)
seq.slab.thick = 0.8*seq.fov(3); % to avoid wrap-around in z
seq.slab.tbw = 12;
seq.slab.type = 'st';  % 'st' = small-tip. 'ex' = 90 degree design
seq.slab.ftype = 'min'; % minimum-phase SLR pulse is well suited for 3D slab excitation
seq.slab.dur = 1.5;  % ms

% SMS excitation pulse parameters
seq.sms.flip = 30;        % flip angle (degrees)
seq.sms.type = 'st';      % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
seq.sms.ftype = 'ls';     
seq.sms.tbw = 6;          % time-bandwidth product
seq.sms.dur = 4;          % msec
seq.sms.sliceSep = seq.fov(3)/mb;   % center-to-center separation between SMS slices (cm)

% slew rates for waveform design (G/cm/ms)
slewRead = [11 15 15];
slewPre = 8;

nx = seq.imSize(1);
ny = seq.imSize(2);
nz = seq.imSize(3);

if mod(ny,seq.Ry) > 0
    error('ny/Ry must be an integer');
end
if mod(nz,seq.Rz) > 0
    error('nz/Rz must be an integer');
end


%% Create modules.txt  and toppeN.entry files

% Entries are tab-separated
fid = fopen('modules.txt', 'wt');
fprintf(fid, 'Total number of unique cores\n');
fprintf(fid, '%d\n', 3);
fprintf(fid, 'fname  duration(us)    hasRF?  hasDAQ?\n');
fprintf(fid, '%s\t0\t1\t0\n', 'tipdown.mod');
fprintf(fid, '%s\t0\t0\t0\n', 'prephase.mod');
fprintf(fid, '%s\t0\t0\t1\n', 'readout.mod');
fclose(fid);

% The entry file can be edited by hand as needed after copying to scanner
fid = fopen('toppeN.entry', 'wt');
fprintf(fid, '/usr/g/research/pulseq/hmri/SMS-EPI/\n');  
fprintf(fid, 'modules.txt\n');
fprintf(fid, 'scanloop.txt\n');
fprintf(fid, '%s\n', 'tipdown.mod');
fprintf(fid, '%s\n', 'readout.mod');
fprintf(fid, 'seqstamp.txt');
fclose(fid);


%% Slab/SMS excitation
if strcmp(scanType, 'SMS')
    tmp = sys.ge;
    tmp.maxSlew = 8;   % G/cm/ms. Reduce PNS during slice rephaser.
    [ex.rf, ex.g, freq] = getsmspulse(seq.sms.flip, seq.voxelSize(3), seq.sms.tbw, seq.sms.dur, ...
        mb, seq.sms.sliceSep, tmp, ...
	    'doSim', true, ...   % Plot simulated SMS slice profile
	    'type', seq.sms.type, ...
	    'ftype', seq.sms.ftype);

    toppe.writemod(sys.ge, 'rf', ex.rf, 'gz', ex.g, ...
	    'ofname', 'tipdown.mod' );

    freq = freq/seq.sms.sliceSep*seq.voxelSize(3); % frequency offset for z shift of seq.voxelSize(3)
else
    % slab select
    toppe.utils.rf.makeslr(seq.slab.flip, seq.slab.thick, ...
        seq.slab.tbw, seq.slab.dur, nz*seq.nSpoilCycles, sys.ge, ...
        'type', seq.slab.type, ...     % 'st' = small-tip. 'ex' = 90 degree design
        'ftype', seq.slab.ftype, ...  
        'spoilDerate', 0.5, ...
        'ofname', 'tipdown.mod');
end


%% Readout 
[gx, gy, gz, gpre, esp, gx1, kz] = getcaipiepireadout(seq.fov, seq.imSize, ...
    seq.Ry, seq.pf_ky, ...
    seq.Rz, seq.Delta, ...
    sys.ge.maxGrad, slewRead, slewPre, ...
    sys.ge.raster*1e3, sys.ge.forbiddenEspRange);

toppe.writemod(sys.ge, ...
    'gx', gpre.x, 'gy', gpre.y, 'gz', gpre.z, ...
    'ofname', 'prephase.mod' );

toppe.writemod(sys.ge, ...
    'gx', gx, 'gy', gy, 'gz', gz, ...
    'ofname', 'readout.mod' );

%% Scan loop
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;

toppe.write2loop('setup', sys.ge, 'version', 4);   % Initialize scanloop.txt

IZ = 1:seq.Rz:nz;

% temporal loop
for ifr = 1:(seq.nCalFrames  + seq.nFrames)

    isCalScan = (ifr < seq.nCalFrames + 1);

    % shift ky sampling pattern (for Ry = 2, this corresponds to UNFOLD)
    %iy = (1 + (-1)^(ifr))/2 + 1;   % 1 or 2
    %a_gy = ((iy-1+0.5)-ny/2)/(ny/2);
    a_gy = (1-isCalScan)*((1-1+0.5)-ny/2)/(ny/2);

    a_gx = (-1)^(ifr+1); % alternate so we can combine frames to simulate flyback EPI

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
        toppe.write2loop('tipdown.mod', sys.ge, ...
            'RFoffset', f, ...
            'RFphase', rfphs);

        % readout 
        % data is stored in 'slice', 'echo', and 'view' indeces
        toppe.write2loop('prephase.mod', sys.ge, ...
            'Gamplitude', [-1.0*a_gx a_gy a_gz]');
        toppe.write2loop('readout.mod', sys.ge, ...
            'Gamplitude', [a_gx (1-isCalScan) (1-isCalScan)]', ...
            'DAQphase', rfphs, ...
            'slice', ii, 'view', ifr);

        % make net gradient area per TR constant (to achieve steady state)
        toppe.write2loop('prephase.mod', sys.ge, ...
            'Gamplitude', [0 -a_gy -a_gz]');

        % update rf phase (RF spoiling)
        rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
        rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
    end
end
toppe.write2loop('finish', sys.ge);
fprintf('\n');

save gx1 gx1
save kz kz

%% Create 'sequence stamp' file for TOPPE.
% This file is listed in the 5th row in toppeN.entry
% NB! The file toppeN.entry must exist in the folder from where this script is called.
toppe.preflightcheck('toppeN.entry', 'seqstamp.txt', sys.ge);

%% create tar file
system('tar cf smsepi.tar toppeN.entry seqstamp.txt scanloop.txt modules.txt *.mod gx1.mat kz.mat');

%tar('epi3d.tar', {entryFile, 'modules.txt', 'scanloop.txt', 'seqstamp.txt', ...
%    'tipdown.mod', 'prephase.mod', 'readout.mod', ...
%    'epi3d.m', 'getparams.m', 'getepireadout.m', ...
%    'gx1.mat', 'kz.mat'});

%fprintf(sprintf('\nPlace %s in /usr/g/research/pulseq/ on scanner host.\n', entryFile));
%fprintf(sprintf('Place all other files in %s on scanner host.\n', filePath));
