function epiInfo = caipiepifmri(scanType, FOV, imSize, Ry, pf_ky, Rz, CaipiShiftZ, flip, nFrames, varargin)
% function epiInfo = caipiepifmri(scanType, FOV, imSize, Ry, pf_ky, Rz, CaipiShiftZ, flip, nFrames, varargin)
%
% Creat 3D/SMS CAIPI EPI fMRI sequence in TOPPE format.
%
% 3D version implements sequence in Narsude et al MRM 2016, 
% "Three-Dimensional Echo Planar Imaging with Controlled Aliasing:
% A Sequence for High Temporal Resolution Functional MRI"
%
% SMS version is similar except partition-encode replaced
% by z shift of excited slices.
%
% Inputs:
%  scanType    string   '3D' or 'SMS'
%  FOV         [1 3]    field of view (cm)
%  imSize      [1 3]    image volume matrix size
%  Ry          [1]      ky acceleration factor (can be non-integer)
%  pf_ky       [1]      Partial Fourier factor (along ky). Range is [0.7 1.0].
%  Rz          [1]      kz acceleration factor (integer)
%  CaipiShiftZ [1]      size of kz step (integer multiples of 1/FOV(3))). Must be > 0.
%  flip        [1 1]    flip angle (degrees)
%  nFrames     [1 1]    number of time frames
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

if strcmp(scanType, 'test')
    sub_test();
    return
end

%% Set defaults for keyword arguments and replace if provided
arg.entryFile = 'toppeN.entry';
arg.filePath = '/usr/g/research/pulseq/fmri/';  % location of scan files on scanner host
arg.rfSpoilSeed = 117;         % RF spoiling phase increment factor (degrees)
arg.mb = Rz;        % sms/multiband factor (number of simultaneous slices)
arg.nCalFrames = 4;   % Odd/even echo calibration frames (ky and kz-encoding off)
arg.doSim = true;     % Plot simulated SMS slice profile
arg.alternateReadout = false;  % flip sign of gx gradient every time-frame
arg.nSpoilCycles = 2;   % number of cycles/voxel of gradient spoiling along z
arg.caipiPythonPath = '~/github/HarmonizedMRI/3DEPI/caipi/';
arg.fatsat       = true;         % add fat saturation pulse?
arg.fatFreqSign = -1;            % sign of fatsat pulse frequency offset
arg.ofname = 'fmri.tar'          % output file name

% EPI train parameters
arg.epiGMax = 5;               % peak gradient amplitude (Gauss/cm)
arg.epiSlewRead = [11 15 15];  % Limit slew rate to this value (Gauss/cm/ms), to control PNS.
arg.epiSlewPre = 10;           % Limit slew rate to this value (Gauss/cm/ms) during prewinder.
arg.forbiddenEspRange = [0.41 0.51];  % forbidden echo spacing range (ms)

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
if imSize(3) > 1  % 3D
    seq.slab.thick = 0.8*FOV(3); % to avoid wrap-around in z
    seq.slab.tbw = 8;
    seq.slab.type = 'st';  % 'st' = small-tip. 'ex' = 90 degree design
    seq.slab.ftype = 'min'; % minimum-phase SLR pulse is well suited for 3D slab excitation
    seq.slab.dur = 1.5;  % ms
else
    seq.slab.thick = FOV(3); 
    seq.slab.tbw = 8;
    seq.slab.type = 'st'; 
    seq.slab.ftype = 'ls'; 
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

seq.voxelSize = FOV./imSize;

nx = imSize(1);
ny = imSize(2);
nz = imSize(3);

if mod(ny, Ry) 
    error('ny/Ry must be an integer');
end
if mod(nz, Rz)
    error('nz/Rz must be an integer');
end


%% Create modules.txt and entryFile

% Entries are tab-separated
nModules = 4 + arg.fatsat;
fid = fopen('modules.txt', 'wt');
fprintf(fid, 'Total number of unique cores\n');
fprintf(fid, '%d\n', nModules);
fprintf(fid, 'fname  duration(us)    hasRF?  hasDAQ?\n');
if arg.fatsat
    fprintf(fid, '%s\t0\t1\t0\n', 'fatsat.mod');
end
fprintf(fid, '%s\t0\t1\t0\n', 'tipdown.mod');
fprintf(fid, '%s\t0\t0\t0\n', 'prephase.mod');
fprintf(fid, '%s\t0\t0\t1\n', 'readout.mod');
fprintf(fid, '%s\t0\t0\t0\n', 'spoiler.mod');
fclose(fid);

% Write entry file.
% This can be edited by hand as needed after copying to scanner.
toppe.writeentryfile(arg.entryFile, ...
    'filePath', arg.filePath);


%% fat sat module
fatsat.flip    = 90;
fatsat.slThick = 1000;       % dummy value (determines slice-select gradient, but we won't use it; just needs to be large to reduce dead time before+after rf pulse)
fatsat.tbw     = 2.0;        % time-bandwidth product
fatsat.dur     = 4.5;        % pulse duration (ms)

b1 = toppe.utils.rf.makeslr(fatsat.flip, fatsat.slThick, fatsat.tbw, fatsat.dur, 1e-6, arg.sys, ...
    'type', 'ex', ...    % fatsat pulse is a 90 so is of type 'ex', not 'st' (small-tip)
    'writeModFile', false);
b1 = toppe.makeGElength(b1);
toppe.writemod(arg.sys, 'rf', b1, 'ofname', 'fatsat.mod', 'desc', 'fat sat pulse');


%% Slab/SMS excitation module
if strcmp(scanType, 'SMS')
    tmp = arg.sys;
    tmp.maxSlew = 8;   % G/cm/ms. Reduce PNS during slice rephaser.
    [ex.rf, ex.g, freq] = getsmspulse(flip, seq.voxelSize(3), seq.sms.tbw, seq.sms.dur, ...
        arg.mb, seq.sms.sliceSep, tmp, ...
	    'doSim', arg.doSim, ...   % Plot simulated SMS slice profile
	    'type', seq.sms.type, ...
	    'ftype', seq.sms.ftype);

    nomflip = 360 * arg.sys.gamma * abs(sum(ex.rf)) * 4e-6  % this seems problematic for SMS pulses

    toppe.writemod(arg.sys, 'rf', ex.rf, 'gz', ex.g, ...
        'nChop', [12 12], ...   % make room for RF ringdown during conversion to Pulseq
        'nomflip', nomflip, ...    % for b1 scaling during conversion to Pulseq
	    'ofname', 'tipdown.mod' );

    freq = freq/seq.sms.sliceSep*seq.voxelSize(3); % frequency offset for z shift of seq.voxelSize(3)
else
    % slab select
    tmp = arg.sys;
    tmp.maxSlew = 10;   % G/cm/ms. Reduce PNS during slice rephaser.
    toppe.utils.rf.makeslr(flip, seq.slab.thick, ...
        seq.slab.tbw, seq.slab.dur, nz*arg.nSpoilCycles, tmp, ...
        'type', seq.slab.type, ...     % 'st' = small-tip. 'ex' = 90 degree design
        'ftype', seq.slab.ftype, ...  
        'spoilDerate', 1.0, ...
        'ofname', 'tipdown.mod');
end


%% EPI readout module
[gx, gy, gz, gpre, esp, gx1, kz, nBlipMax] = getcaipiepireadout(FOV, imSize, ...
    Ry, pf_ky, Rz, CaipiShiftZ, ...
    'gMax', arg.epiGMax, ...
    'slewRead', arg.epiSlewRead, ...
    'slewPre', arg.epiSlewPre, ...
    'caipiPythonPath', arg.caipiPythonPath, ...
    'fbesp', arg.forbiddenEspRange);

toppe.writemod(arg.sys, ...
    'gx', gpre.x, 'gy', gpre.y, 'gz', gpre.z, ...
    'ofname', 'prephase.mod' );

toppe.writemod(arg.sys, ...
    'gx', gx, 'gy', gy, 'gz', gz, ...
    'ofname', 'readout.mod' );


%% Spoiler (gradient crusher) module
gspoil = toppe.utils.makecrusher(arg.nSpoilCycles, seq.voxelSize(1), arg.sys, 0, ...
    0.7*arg.sys.maxSlew/sqrt(2), arg.sys.maxGrad/sqrt(2));
toppe.writemod(arg.sys, 'ofname', 'spoiler.mod', ...
    'gx', gspoil, 'gz', gspoil);


%% Scan loop
rfphs = 0;              % radians
rfSpoilSeed_cnt = 0;

toppe.write2loop('setup', arg.sys, 'version', 4);   % Initialize scanloop.txt

% interleaved slice ordering for SMS
IZ = 1:Rz:nz; 
if strcmp(scanType, 'SMS')
    Nex = imSize(3)/Rz;  % number of excitations to cover one image volume
    if ~mod(Nex,2)
        warning('For better interleaving performance, imSize(3)/Rz should be odd');
    end
    II = [1:2:Nex 2:2:Nex];
else
    II = 1:length(IZ);
end

% temporal loop
for ifr = 1:(arg.nCalFrames  + nFrames)
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bFrame %d of %d', ifr, arg.nCalFrames + nFrames);
    isCalScan = (ifr < arg.nCalFrames + 1);

    % y readout prephaser gradient amplitude scaling (normalized)
    a_gy = (1-isCalScan)*((1-1+0.5)-ny/2)/(ny/2);

    if ~isCalScan & arg.alternateReadout
        a_gx = (-1)^(ifr+1); % alternate so we can combine frames to simulate flyback EPI
    else
        a_gx = 1;
    end

    % z encoding / SMS slice shift loop
    for ii = II

        if strcmp(scanType, 'SMS')
		    f = round((ii-0.5-length(IZ)/2)*freq);  % frequency offset (Hz) for slice shift
            a_gz = 0;  % prephaser amplitude
        else
            f = 0;
            a_gz = (1-isCalScan)*((IZ(ii)-1+0.5)-nz/2)/(nz/2);
        end

        % fat sat
        if(arg.fatsat)
            fatChemShift = 3.5;  % fat/water chemical shift (ppm)
            fatFreq = arg.fatFreqSign*arg.sys.gamma*1e4*arg.sys.B0*fatChemShift*1e-6;  % Hz
            toppe.write2loop('fatsat.mod', arg.sys, ...
                'RFoffset', round(fatFreq), ...   % Hz
                'RFphase', rfphs);         % radians

            toppe.write2loop('spoiler.mod', arg.sys);

            rfphs = rfphs + (arg.rfSpoilSeed/180*pi)*rfSpoilSeed_cnt ;  % radians
            rfSpoilSeed_cnt = rfSpoilSeed_cnt + 1;
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

        % spoiler gradient
        toppe.write2loop('spoiler.mod', arg.sys);

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

%% Return struct needed for recon
epiInfo.imSize   = imSize;
epiInfo.pf_ky    = pf_ky;
epiInfo.Ry       = Ry;
epiInfo.FOV      = FOV;
epiInfo.gpre     = gpre;  % x prewinder (G/cm)
epiInfo.gx       = gx;      % full echo train (G/cm)
epiInfo.gx1      = gx1;    % one echo
epiInfo.esp      = esp;
epiInfo.nBlipMax = nBlipMax;  % number of samples durings turns
save epiInfo epiInfo

%% create tar file
system(sprintf('tar cf %s %s seqstamp.txt scanloop.txt modules.txt *.mod epiInfo.mat', arg.ofname, arg.entryFile));

toppe.utils.scanmsg(arg.entryFile);


return




%% test function / example usage
function sub_test()

    % System hardware specs
    sysGE = toppe.systemspecs('maxSlew', 20, 'slewUnit', 'Gauss/cm/ms', ...
        'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
        'myrfdel', 200, ...  % psd_rf_wait
        'daqdel', 100, ...   % psd_grd_wait
        'timessi', 100, ...
        'gradient', 'hrmb');  

    % sequence parameters
    res = 0.275;  % cm
    imSize = [80 80 60];
    FOV = imSize*res;  % cm
    Ry = 1;
    pf_ky = 0.8;
    Rz = 4;
    CaipiShiftZ = 2;
    flip = 15;
    nFrames = 10;
    
    % create 3D EPI fMRI sequence. Writes TOPPE files to the current folder.
    caipiepifmri('3D', FOV, imSize, Ry, pf_ky, Rz, CaipiShiftZ, flip, nFrames, ...
        'sys', sysGE, ...
        'fatsat', true, ...
        'epiGMax', 5, ...
        'epiSlewRead', [11 15 15], ...
        'epiSlewPre', 10, ...
        'forbiddenEspRange', [0.41 0.51]);

    % Display scan loop
    % The 'playseq' function acts a lot like the interpreter, i.e.,
    % it reads modules.txt and scanloop.txt in the current working directory
    % and 'executes' the sequence.
    % The displayed sequence timing (determined by sys) tries to exactly match
    % what you'd see on the scanner.
    nModsPerTR = 5;
    reply = input('Display 3D EPI scan loop? Y/N [Y] ', 's');
    if isempty(reply) | strcmp(upper(reply), 'Y')
        toppe.playseq(nModsPerTR, sysGE, 'nTRskip', 1, 'tpause', 0.1);
    end
    
    % Create SMS EPI fMRI sequence
    caipiepifmri('SMS', FOV, imSize, Ry, pf_ky, Rz, CaipiShiftZ, flip, nFrames, ...
        'sys', sysGE, ...
        'fatsat', true, ...
        'epiGMax', 5, ...
        'epiSlewRead', [11 15 15], ...
        'epiSlewPre', 10, ...
        'forbiddenEspRange', [0.41 0.51]);

    % Display scan loop
    nModsPerTR = 5;
    reply = input('Display SMS EPI scan loop? Y/N [Y] ', 's');
    if isempty(reply) | strcmp(upper(reply), 'Y')
        toppe.playseq(nModsPerTR, sysGE, 'nTRskip', 1, 'tpause', 0.1);
    end

return
