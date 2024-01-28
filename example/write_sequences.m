% create SMS EPI sequence files, in the following order:
% 1. mb=6 sequence. Defines readout and adc events
% 2. calibration scan (blips off, for receive gain calibration and EPI ghost correction)
% 3. mb=1 single-slice 'ACS' reference scan

TODO = [1 1 1];

sysGE = toppe.systemspecs();  % for plotting

fatSat = false;
RFspoil = false;

% mb=6 sequence. Defines readout and adc events
% Design sub-sequence containing 40 shots = 4 frames 
% We choose 4 shots since RF spoiling phase (for rf_inc = 117) repeats every 80 RF shots
% (fat sat also spoils so only need 40 TRs not 80)
% But RF spoiling anyhow probably isn't doing much since TR>>T2.
mb = 6; Ry = 1; Rz = mb; caipiShiftZ = 2;
nDummyFrames = 0;
nFrames = 1;
if TODO(1)
    [gro, adc] = writeEPI(voxelSize, [nx ny nz], TE, TR, alpha, mb, pf_ky, Ry, Rz, caipiShiftZ, nFrames, nDummyFrames, 'SMS', ...
        'seqName', 'mb6', ...
        'fatSat', fatSat, ...
        'RFspoil', RFspoil, ...
        'doRefScan', false, ...
        'simulateSliceProfile', true);
    dur = 5*40 + 10;    % scan duration, sec
    nRuns = dur/(nFrames*TR)     % CV8/opuser8 on UI
end

% Prescan (Rx gain) and EPI ghost calibration
% (reuse mb, Ry, etc from above)
nDummyFrames = 2;   % to reach steady state
nFrames = 2;        % 1st frame is ghost calibration (y blips off)
if TODO(2)
    writeEPI(voxelSize, [nx ny nz], TE, TR, alpha, mb, pf_ky, Ry, Rz, caipiShiftZ, nFrames, nDummyFrames, 'SMS', ...
        'seqName', 'cal', ...
        'fatSat', fatSat, ...
        'RFspoil', RFspoil, ...
        'doRefScan', true, ...
        'gro', gro, 'adc', adc);
end

% 2D mb=1 sequence for reference, and for distortion-matched sensitivity maps
mb = 1; Ry = 1; Rz = mb; caipiShiftZ = 0;
nDummyFrames = 0;
nFrames = 1; 
pf_ky_ref = 1.0;
if TODO(3)
    writeEPI(voxelSize, [nx ny nz], TE, 5*TR, alpha, mb, pf_ky_ref, Ry, Rz, caipiShiftZ, nFrames, nDummyFrames, 'SMS', ...
        'seqName', '2d', ...
        'fatSat', fatSat, ...
        'RFspoil', RFspoil, ...
        'doRefScan', false, ...
        'simulateSliceProfile', true, ...
        'gro', gro, 'adc', adc);
end

% 3D mb=1 sequence for reference, and for distortion-matched sensitivity maps
mb = 1; Ry = 1; Rz = mb; caipiShiftZ = 0;
nDummyFrames = 1;
nFrames = 2;      % first frame is for Rx gain calibration (APS)
pf_ky_ref = 1.0;
if TODO(4)
    writeEPI(voxelSize, [nx ny nz], TE, 5*TR, alpha, 1, pf_ky_ref, Ry, Rz, caipiShiftZ, nFrames, nDummyFrames, '3D', ...
        'seqName', '3d', ...
        'fatSat', fatSat, ...
        'RFspoil', RFspoil, ...
        'doRefScan', true, ...
        'gro', gro, 'adc', adc);
end

% 3D GRE for B0 and sensitivity maps
if TODO(5)
    writeB0;
end

system('tar cf all.tar *.tar *.seq *.m');

return


