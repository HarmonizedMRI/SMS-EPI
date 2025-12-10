% Create multi-echo SMS-EPI sequence files.
% 1. ME SMS-EPI fMRI sequence
% 2. noise scan (same as #1 except no rf)
% 3. EPI ghosting and receive gain calibration scan (blips off)
% 4. mb=1, Ry=2, ME slice grappa calibration scan
% 5. mb=1, Ry=1 grappa calibration scan
% 6. 3D GRE B0 mapping (also usable for sensitivity maps or grappa calibration)

system('rm -f *.tar *.seq *.pge *.entry');

TODO = ones(1,6);

%----------------------------------------------------------
% sequence parameters
%----------------------------------------------------------
% SMS=4, Ry=2, pf_ky = 0.85, fov 24x24cm, 80x80, 3mm iso:
% Cohen, Alexander D., et al. "Detecting task functional MRI activation 
% using the multiband multiecho (MBME) echo‐planar imaging (EPI) sequence." 
% Journal of Magnetic Resonance Imaging 53.5 (2021): 1366-1374.

fatSat = true;
RFspoil = true;
nx = 80; ny = nx; nz = 44;
voxelSize = [3 3 3] * 1e-3;   % m
alpha = 55;
pf_ky = 0.85; %(nx-3*6)/nx;

TR = 0.9;  %0.1*11;   % volume TR (sec)

PNSwt = [0.8 1 0.7];   % PNS direction/channel weights

mb = 4; Ry = 2; caipiShiftZ = 2;
[IY, IZ] = getcaipi(ny, nz, Ry, mb, caipiShiftZ, '3DEPI/caipi');
etl = 2*ceil(pf_ky*ny/Ry/2);  % echo train length for one TE 
IY = IY(end-etl+1:end);
IZ = IZ(end-etl+1:end);
nTE = 3;   
IY = repmat(IY, nTE, 1);
IZ = repmat(IZ, nTE, 1);

%----------------------------------------------------------
% scanner hardware settings
%----------------------------------------------------------
vendor = 'GE';   % 'GE' or 'Siemens'

switch lower(vendor)
    case 'ge'
        % Settings for GE scanners at U-Mich
        gySign = 1;
        freqSign = 1;
        fatFreqSign = -1;
        doConj = false;
        gradRasterTime = 4e-6;
        rfRasterTime = 4e-6;
        blockDurationRaster = 4e-6;
        B0 = 3;
        segmentRingdownTime = 117e-6;

        psd_rf_wait = 100e-6;  % RF-gradient delay, scanner specific (s)
        psd_grd_wait = 100e-6; % ADC-gradient delay, scanner specific (s)
        b1_max = 0.25;         % Gauss
        g_max = 5;             % Gauss/cm
        slew_max = 20;         % Gauss/cm/ms
        coil = 'xrm';          % 'hrmbuhp' (UHP); 'xrm' (MR750)
        sysGE = pge2.opts(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, coil);

        % settings for writing pge<>.entry files
        pgeFilePath = '/srv/nfs/psd/usr/psd/pulseq/v7/sequences/hfmri';
        entryFileNumber = 410;  % starting value

    case 'siemens'
        % Settings for Siemens Vida bay MR7 @ U-Mich University Hospital
        gySign = 1;
        freqSign = -1;
        fatFreqSign = -1;  % confirmed on Prisma Fit at VA Tech 2024-Oct-24, and Urbana-Champaign
        doConj = false;
        gradRasterTime = 10e-6;
        rfRasterTime = 1e-6;
        blockDurationRaster = 10e-6;
        B0 = 2.89;
        segmentRingdownTime = 0;
end

% slew is set to avoid forbidden EPI echo spacings for MR750 and UHP scanners
sys = mr.opts('maxGrad', 40, 'gradUnit','mT/m', ...
              'maxSlew', 120, 'slewUnit', 'T/m/s', ...
              'rfDeadTime', 100e-6, ...
              'rfRingdownTime', 60e-6, ...
              'adcDeadTime', 0e-6, ...
              'adcRasterTime', 4e-6, ...
              'gradRasterTime', gradRasterTime, ...
              'rfRasterTime', rfRasterTime, ...
              'blockDurationRaster', blockDurationRaster, ...
              'B0', B0);

%----------------------------------------------------------
% options passed to writeEPI.m
%----------------------------------------------------------
opts = struct('fatSat', fatSat, ...
    'RFspoil', RFspoil, ...
    'gySign', gySign, ...
    'freqSign', freqSign, ...
    'fatFreqSign', fatFreqSign, ...
    'doConj', doConj, ...
    'segmentRingdownTime', segmentRingdownTime, ...
    'doRefScan', false, ...
    'trigOut', false, ...
    'doNoiseScan', false, ...
    'plot', false, ...
    'simulateSliceProfile', false, ...
    'echo', []);

%----------------------------------------------------------
% write sequences
%----------------------------------------------------------

nFrames = 1;

% fmri.seq 
% Determines readout gradient and ADC event (defined in return struct echo)
if TODO(1)
    fn = 'fmri';
    nFramesTmp = 4;   % since opnex is limited
    [~, echo] = writeEPI(fn, sys, voxelSize, [nx ny nz], TR, alpha, mb, IY, IZ, ...
        nFramesTmp, 'SMS', opts);
    if strcmp(lower(vendor), 'ge')
        pge2.seq2ge(fn, sysGE, length(IY)*nz/mb, PNSwt);
        pge2.writeentryfile(entryFileNumber, fn, 'path', pgeFilePath);
    end
end

% noise.seq
if TODO(2)
    fn = 'noise';
    writeEPI(fn, sys, voxelSize, [nx ny nz], TR, alpha, mb, IY, IZ, ...
        nFrames, 'SMS', ...
        pge2.utils.setfields(opts, ...
            'doNoiseScan', true,  ...
            'echo', echo) );
    if strcmp(lower(vendor), 'ge')
        pge2.seq2ge(fn, sysGE, 1, PNSwt);
        entryFileNumber = entryFileNumber + 1;
        pge2.writeentryfile(entryFileNumber, fn, 'path', pgeFilePath);
    end
end

% epical.seq
%  - EPI ghost calibration
%  - Receive gain calibration (GE)
if TODO(3)
    fn = 'epical';
    writeEPI(fn, sys, voxelSize, [nx ny nz], 2*TR*mb, alpha, 1, IY, IZ, ...
        nFrames, 'SMS', ...
        pge2.utils.setfields(opts, 'doRefScan', true, 'echo', echo));
    if strcmp(lower(vendor), 'ge')
        entryFileNumber = entryFileNumber + 1;
        pge2.writeentryfile(entryFileNumber, fn, 'path', pgeFilePath);
        pislquant = round(length(IY)*nz/2);
        pge2.seq2ge(fn, sysGE, pislquant, PNSwt);
    end
end

% slgcal.seq
%  - slice GRAPPA calibration
if TODO(4)
    fn = 'slgcal';
    writeEPI(fn, sys, voxelSize, [nx ny nz], 2*TR*mb, alpha, 1, IY, ones(size(IZ)), ...
        nFrames, 'SMS', ...
        pge2.utils.setfields(opts, 'echo', echo));
    if strcmp(lower(vendor), 'ge')
        pge2.seq2ge(fn, sysGE, 1, PNSwt);
        entryFileNumber = entryFileNumber + 1;
        pge2.writeentryfile(entryFileNumber, fn, 'path', pgeFilePath);
    end
end

% grappacal.seq
%  - GRAPPA calibration
if TODO(5)
    fn = 'grappacal';
    [IYtmp, IZtmp] = getcaipi(ny, nz, 1, 1, 0, '3DEPI/caipi');
    writeEPI(fn, sys, voxelSize, [nx ny nz], 2*TR*mb, alpha, 1, IYtmp, IZtmp, ...
        nFrames, 'SMS', ...
        pge2.utils.setfields(opts, 'echo', echo));
    if strcmp(lower(vendor), 'ge')
        pge2.seq2ge(fn, sysGE, 1, PNSwt);
        entryFileNumber = entryFileNumber + 1;
        pge2.writeentryfile(entryFileNumber, fn, 'path', pgeFilePath);
    end
end

% b0.seq
%  - 3D B0 field / coil sensitivity maps
if TODO(6)
    fn = 'b0';
    writeB0(fn, ...
        pge2.utils.setfields(sys, ...
            'rfRingdownTime', 300e-6, ... % make room for psd_rf_wait (GE)    
            'maxSlew', sys.maxSlew*0.5, ...
            'maxGrad', sys.maxGrad*0.5), ...
        voxelSize, [nx nx nx], 4);
    if strcmp(lower(vendor), 'ge')
        pge2.seq2ge(fn, sysGE, ny, PNSwt);
        entryFileNumber = entryFileNumber + 1;
        pge2.writeentryfile(entryFileNumber, fn, 'path', pgeFilePath);
    end
end

if strcmp(lower(vendor), 'ge')
    ofn = 'pgescans-' + string(date) + '.tar';
    system(sprintf('rm %s', ofn));
    system(sprintf('tar cf %s *.entry *.pge', ofn));
end
    
