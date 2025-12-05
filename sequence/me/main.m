% Create multi-echo SMS-EPI sequence files.
% 1. ME SMS-EPI fMRI sequence. Defines readout (gro) and adc events
% 2. EPI ghosting and receive gain calibration scan (blips off)
% 3. mb=1, Ry=2, ME slice grappa calibration scan
% 4. mb=1, Ry=1 grappa calibration scan
% 5. 3D GRE B0 mapping (can also be used for sens maps)

TODO = [1 0 0 0 0];

vendor = 'GE';   % 'GE' or 'Siemens'
fatSat = true;
RFspoil = true;

% Acquisition parameters
% e.g., SMS=4, Ry=2, pf_ky = 0.85, fov 24x24cm, 80x80, 3mm iso:
% Cohen, Alexander D., et al. "Detecting task functional MRI activation 
% using the multiband multiecho (MBME) echo‐planar imaging (EPI) sequence." 
% Journal of Magnetic Resonance Imaging 53.5 (2021): 1366-1374.
nx = 80; ny = nx; nz = 44;
voxelSize = [3 3 3] * 1e-3;   % m
alpha = 60;
pf_ky = 0.85; %(nx-3*6)/nx;

TR = 0.1*11;   % volume TR (sec)

switch lower(vendor)
    case 'ge'
        % Settings for GE scanners at U-Mich
        gySign = 1;
        freqSign = 1;
        fatFreqSign = -1;
        doConj = false;
    case 'siemens'
        % Settings for Siemens Vida bay MR7 @ U-Mich University Hospital
        gySign = 1;
        freqSign = -1;
        fatFreqSign = -1;  % Yes, confirmed on Prisma Fit at VA Tech 2024-Oct-24, and Urbana-Champaign
        doConj = false;
end

% input options to writeEPI.m
opts = struct('vendor', vendor, ...
    'seqName', sprintf('mb%d', mb), ...
    'fatSat', fatSat, ...
    'RFspoil', RFspoil, ...
    'gySign', gySign, ...
    'freqSign', freqSign, ...
    'fatFreqSign', fatFreqSign, ...
    'doConj', doConj, ...
    'doRefScan', false, ...
    'trigOut', false, ...
    'doNoiseScan', false, ...
    'simulateSliceProfile', true);

% SMS-EPI sequence
mb = 4; Ry = 2; caipiShiftZ = 2;
[IY, IZ] = getcaipi(ny, nz, Ry, mb, caipiShiftZ, '3DEPI/caipi');
etl = 2*ceil(pf_ky*ny/Ry/2);  % echo train length. even
IY = IY(end-etl+1:end);
IZ = IZ(end-etl+1:end);
nTE = 3;   
IY = repmat(IY, nTE, 1);
IZ = repmat(IZ, nTE, 1);
nFrames = 1;
if TODO(1)
    writeEPI('fmri', voxelSize, [nx ny nz], TR, alpha, mb, IY, IZ, nFrames, 'SMS', opts);

    dur = 5*40 + 10;    % scan duration, sec
    dur = 5*60 + 10;    % scan duration, sec
    nRuns = dur/(nFrames*TR);     % CV8/opuser8 on UI
end

% For prescan (Rx gain) and EPI ghost calibration
% (reuse mb, Ry, etc from above)
if TODO(2)
    opts.doRefScan = true;
    writeEPI('cal', voxelSize, [nx ny nz], TR, alpha, mb, IY, IZ, nFrames, 'SMS', opts);
    opts.doRefScan = false;
end

% 2D mb=1 sequence for:
%  slice GRAPPA reference 
%  slice-by-slice EPI calibration
%  having receiving being set automatically during prescan (GE)
mb = 1; Ry = 1; Rz = mb; caipiShiftZ = 0;
nDummyFrames = 0;
nFrames = 2; 
if TODO(3)
    writeEPI('2d', voxelSize, [nx ny nz], TE, 2*6*TR, alpha, mb, pf_ky, Ry, Rz, caipiShiftZ, nFrames, nDummyFrames, 'SMS', opts);
end

% noise scan
if TODO(4)
    opts.doNoiseScane = true;
    writeEPI('noise', voxelSize, [nx ny nz], TE, TR, alpha, mb, pf_ky, Ry, Rz, caipiShiftZ, nFrames, nDummyFrames, 'SMS', opts);
    opts.doNoiseScane = false;
end

% 3D GRE for B0 and sensitivity maps
if TODO(5)
    writeB0;
end
