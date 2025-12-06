% Create multi-echo SMS-EPI sequence files.
% 1. ME SMS-EPI fMRI sequence. Defines readout (gro) and adc events
% 2. noise scan (same as 1. except no rf)
% 3. EPI ghosting and receive gain calibration scan (blips off)
% 4. mb=1, Ry=2, ME slice grappa calibration scan
% 5. mb=1, Ry=1 grappa calibration scan
% 6. 3D GRE B0 mapping (can also be used for sens maps and/or grappa calibration)

TODO = [1 1 1 1 1 0];
TODO = [1 1 0 0 0 1];
TODO = ones(1,6);


%----------------------------------------------------------
% Sequence parameters
%----------------------------------------------------------
% e.g., SMS=4, Ry=2, pf_ky = 0.85, fov 24x24cm, 80x80, 3mm iso:
% Cohen, Alexander D., et al. "Detecting task functional MRI activation 
% using the multiband multiecho (MBME) echo‐planar imaging (EPI) sequence." 
% Journal of Magnetic Resonance Imaging 53.5 (2021): 1366-1374.

fatSat = true;
RFspoil = true;
nx = 80; ny = nx; nz = 44;
voxelSize = [3 3 3] * 1e-3;   % m
alpha = 60;
pf_ky = 0.85; %(nx-3*6)/nx;

TR = 0.1*11;   % volume TR (sec)

mb = 4; Ry = 2; caipiShiftZ = 2;
[IY, IZ] = getcaipi(ny, nz, Ry, mb, caipiShiftZ, '3DEPI/caipi');
etl = 2*ceil(pf_ky*ny/Ry/2);  % echo train length. even
IY = IY(end-etl+1:end);
IZ = IZ(end-etl+1:end);
nTE = 3;   
IY = repmat(IY, nTE, 1);
IZ = repmat(IZ, nTE, 1);


%----------------------------------------------------------
% scanner settings
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
    case 'siemens'
        % Settings for Siemens Vida bay MR7 @ U-Mich University Hospital
        gySign = 1;
        freqSign = -1;
        fatFreqSign = -1;  % Yes, confirmed on Prisma Fit at VA Tech 2024-Oct-24, and Urbana-Champaign
        doConj = false;
        gradRasterTime = 10e-6;
        rfRasterTime = 4e-6;
        blockDurationRaster = 10e-6;
        B0 = 2.89;
end

sys = mr.opts('maxGrad', 40, 'gradUnit','mT/m', ...
              'maxSlew', 150, 'slewUnit', 'T/m/s', ...
              'rfDeadTime', 100e-6, ...
              'rfRingdownTime', 60e-6, ...
              'adcDeadTime', 0e-6, ...
              'adcRasterTime', 4e-6, ...
              'gradRasterTime', gradRasterTime, ...
              'rfRasterTime', rfRasterTime, ...
              'blockDurationRaster', blockDurationRaster, ...
              'B0', B0);

%----------------------------------------------------------
% input options for writeEPI.m
%----------------------------------------------------------
opts = struct('fatSat', fatSat, ...
    'RFspoil', RFspoil, ...
    'gySign', gySign, ...
    'freqSign', freqSign, ...
    'fatFreqSign', fatFreqSign, ...
    'doConj', doConj, ...
    'doRefScan', false, ...
    'trigOut', false, ...
    'doNoiseScan', false, ...
    'plot', true, ...
    'simulateSliceProfile', true);

%----------------------------------------------------------
% write sequences
%----------------------------------------------------------

% fmri.seq 
nFrames = 1;
if TODO(1)
    writeEPI('fmri', sys, voxelSize, [nx ny nz], TR, alpha, mb, IY, IZ, nFrames, 'SMS', opts);
end

% noise.seq
if TODO(2)
    writeEPI('noise', sys, voxelSize, [nx ny nz], TR, alpha, mb, IY, IZ, nFrames, 'SMS', ...
        pge2.utils.override(opts, 'doNoiseScan', true));
end

% epical.seq
%  - EPI ghost calibration
if TODO(3)
    writeEPI('epical', sys, voxelSize, [nx ny nz], TR*mb, alpha, 1, IY, IZ, nFrames, 'SMS', ...
        pge2.utils.override(opts, 'doRefScan', true));
end

% slgcal.seq
%  - slice GRAPPA calibration
if TODO(4)
    writeEPI('slgcal', sys, voxelSize, [nx ny nz], TR*mb, alpha, 1, IY, ones(size(IZ)), nFrames, 'SMS', opts);
end

% grappacal.seq
%  - GRAPPA calibration
if TODO(5)
    [IYtmp, IZtmp] = getcaipi(ny, nz, 1, 1, 0, '3DEPI/caipi');
    writeEPI('grappacal', sys, voxelSize, [nx ny nz], TR*mb, alpha, 1, IYtmp, IZtmp, nFrames, 'SMS', opts);
end

% b0.seq
%  - 3D B0 field / coil sensitivity maps
if TODO(6)
    writeB0('b0', ...
        pge2.utils.override(sys, 'maxSlew', sys.maxSlew*0.5, 'maxGrad', sys.maxGrad*0.5), ...
        voxelSize, [nx nx nx], 4);
end

