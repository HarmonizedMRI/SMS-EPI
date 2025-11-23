function [gro, adc] = writeEPI(voxelSize, N, TE, TR, alpha, mb, pf_ky, Ry, Rz, caipiShiftZ, nFrames, nDummyFrames, tp, varargin)
% function [gro, adc] = writeEPI(voxelSize, N, TE, TR, alpha, mb, pf_ky, Ry, Rz, caipiShiftZ, nFrames, nDummyFrames, tp, varargin)
%
% SMS-EPI sequence in Pulseq
%
% Inputs:
%   voxelSize    [1 3]    meters
%   N            [1 3]    image matrix size 
%   TE           [1]      sec
%   alpha        [1]      flip angle (degrees)
%   mb           [1]      multiband/SMS factor
%   pf_ky        [1]      partial Fourier factor
%   Ry           [1]      ky undersampling factor
%   Rz           [1]      kz undersampling factor
%   caipiShiftZ  [1]      caipi shift along z for every ky blip
%   nFrames      [1]      number of time frames (image volumes)
%   nDummyFrames [1]      number of frames w/o data acquisition to reach steady state
%   tp         string   'SMS', 'SE', or '3D'. If 'SE', mb is set to 1
%
% Keyword-argument input options:
%   gro          struct   readout gradient struct. Default: create a new one and return from this function
%   adc          struct   acquisition struct. Default: create a new one and return from this function
%
% Outputs:
%   gro 
%   adc

tp = 'ME';

[nx ny nz] = deal(N(1), N(2), N(3));

np = nz/mb;   % number of excitations/partitions (sets of SMS slices)

fprintf('mb=%d\n', mb); 

% parse input options
arg.gro = [];
arg.adc = [];
arg.fatSat = true;  
arg.RFspoil = true;
arg.simulateSliceProfile = false;    % simulate SMS profile and display
arg.gzPreOn = true;                  % prephase along kz by -mb/2*deltak
arg.plot = false;
arg.segmentRingdownTime = 120e-6;    % segment ringdown time for Pulseq on GE 
arg.doRefScan = false;               % do EPI ghost reference scan in frame 1
arg.doNoiseScan = false;
arg.toGE = true;
arg.seqName = 'epi';
arg.gySign = +1;
arg.freqSign = +1;
arg.fatFreqSign = -1;
arg.doConj = false;
arg.trigOut = false;
%arg.caipiPythonPath = '~/github/HarmonizedMRI/3DEPI/caipi/';
arg.caipiPythonPath = '3DEPI/caipi/';

arg = toppe.utils.vararg_pair(arg, varargin);

if arg.doNoiseScan
    arg.fatSat = false;
end

if mb == 1
%    arg.gzPreOn = false;
end

warning('OFF', 'mr:restoreShape');

%% Define experimental parameters
sys = mr.opts('maxGrad', 40, 'gradUnit','mT/m', ...
              'maxSlew', 150, 'slewUnit', 'T/m/s', ...
              'rfDeadTime', 0e-6, ...
              'rfRingdownTime', 0e-6, ...
              'adcDeadTime', 0e-6, ...
              'adcRasterTime', 4e-6, ...
              'gradRasterTime', 4e-6, ...
              'rfRasterTime', 4e-6, ...
              'blockDurationRaster', 4e-6, ...
              'B0', 3.0);

% reduced slew for spoilers
sys2 = mr.opts('maxGrad', 40, 'gradUnit','mT/m', ...
              'maxSlew', 60, 'slewUnit', 'T/m/s', ...
              'rfDeadTime', sys.rfDeadTime, ...
              'rfRingdownTime', sys.rfRingdownTime, ...
              'adcDeadTime', sys.adcDeadTime, ...
              'gradRasterTime', sys.gradRasterTime, ...
              'rfRasterTime', sys.rfRasterTime, ...
              'blockDurationRaster', sys.blockDurationRaster, ...
              'B0', sys.B0);

fov = voxelSize .* [nx ny nz];       % FOV (m)

if strcmp(tp, '3D')
    slThick = 0.85*fov(3);
else
    slThick = fov(3)/nz;
end

dwell = 4e-6;                    % ADC sample time (s)

etl = 2*ceil(pf_ky*ny/Ry/2);   % echo train length. even

nCyclesSpoil = 2;    % number of spoiler cycles, along x and z
rfSpoilingInc = 117;                % RF spoiling increment (degrees)

if strcmp(tp, '3D')
    rfTB  = 8;          % RF pulse time-bandwidth product
    rfDur = 4e-3;       % RF pulse duration (s)
else
    rfTB  = 6;
    rfDur = 8e-3;
end

fatChemShift = 3.5*1e-6;                        % 3.5 ppm
fatOffresFreq = sys.gamma*sys.B0*fatChemShift;  % Hz


%% Create fat sat pulse 

fatsat.flip    = 90;      % degrees
fatsat.slThick = 1e5;     % dummy value (determines slice-select gradient, but we won't use it; just needs to be large to reduce dead time before+after rf pulse)
fatsat.tbw     = 3.5;     % time-bandwidth product
fatsat.dur     = 8.0;     % pulse duration (ms)

% RF waveform in Gauss
wav = toppe.utils.rf.makeslr(fatsat.flip, fatsat.slThick, fatsat.tbw, fatsat.dur, 1e-6, toppe.systemspecs(), ...
    'tp', 'ex', ...    % fatsat pulse is a 90 so is of type 'ex', not 'st' (small-tip)
    'ftype', 'min', ...
    'writeModFile', false);

% Convert from Gauss to Hz, and interpolate to sys.rfRasterTime
rfp = rf2pulseq(wav, 4e-6, sys.rfRasterTime);

% Create pulseq object
% Try to account for the fact that makeArbitraryRf scales the pulse as follows:
% signal = signal./abs(sum(signal.*opt.dwell))*flip/(2*pi);
flip = fatsat.flip/180*pi;
flipAssumed = abs(sum(rfp));
rfsat = mr.makeArbitraryRf(rfp, ...
    flip*abs(sum(rfp*sys.rfRasterTime))*(2*pi), ...
    'use', 'excitation', ...
    'system', sys);
rfsat.signal = rfsat.signal/max(abs(rfsat.signal))*max(abs(rfp)); % ensure correct amplitude (Hz)
rfsat.freqOffset = arg.fatFreqSign*425;  % Hz


%% excitation pulse
sysGE = toppe.systemspecs('maxGrad', sys.maxGrad/sys.gamma*100, ...   % G/cm
    'maxSlew', 0.7*sys.maxSlew/sys.gamma/10, ...           % G/cm/ms
    'maxRF', 0.15);
sliceSep = fov(3)/mb;   % center-to-center separation between SMS slices (m)
% freq = frequency offset (Hz) corresponding to a sliceSep
if strcmp(tp, 'SMS') | strcmp(tp, 'SE')
    ftype = 'ls';
else
    ftype = 'min';
end
[rf, gzRF, freq, t_rf_center] = getsmspulse(alpha, slThick, rfTB, rfDur, ...
    mb, sliceSep, sysGE, sys, ...
    'doSim', arg.simulateSliceProfile, ...    % Plot simulated SMS slice profile
    'type', 'st', ...     % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
    'noRfOffset', strcmp(tp, '3D'), ...   % don't shift slice (slab) for 3D
    'ftype', ftype);      % filter design. 'ls' = least squares

rf.use = 'excitation';

freq = arg.freqSign * freq;

if arg.doConj
    rf.signal = conj(rf.signal);
end


%% Get CAIPI sampling pattern (for one shot/echo train)

% create caipi.mat, and load it
pyFile = [arg.caipiPythonPath 'skippedcaipi_sampling.py'];
pyCmd = sprintf('python3 %s %d %d %d %d %d %d', ...
    pyFile, ny, nz, Ry, Rz, caipiShiftZ, 1);

if system(pyCmd) ~= 0
    fprintf('Open a terminal and run the following python command:\n\t%s\n', pyCmd);
    input('\tWhen done, press Enter to continue');
end

load caipi

% kz and ky indeces (multiples of deltak)
kyInds = double(indices((end-etl+1):end, 2));
kzInds = double(indices((end-etl+1):end, 1));

% ky/kz encoding blip amplitude along echo train (multiples of deltak)
kyStep = diff(kyInds);
kzStep = diff(kzInds);

kyStepMax = max(abs(kyStep));
kzStepMax = max(abs(kzStep));


%% Define readout gradients and ADC event
% The Pulseq toolbox really shines here!

deltak = 1./fov;

% Start with the blips
commonRasterTime = 20e-6;   
gyBlip = mr.makeTrapezoid('y', sys, 'Area', kyStepMax*deltak(2));
gzBlip = mr.makeTrapezoid('z', sys, 'Area', kzStepMax*deltak(3));

blipDuration = max(mr.calcDuration(gyBlip), mr.calcDuration(gzBlip));
maxBlipArea = max(gyBlip.area, gzBlip.area);

% Readout trapezoid
if isempty(arg.gro) 
    systmp = sys;
    systmp.maxGrad = deltak(1)/dwell;  % to ensure >= Nyquist sampling
    gro = trap4ge(mr.makeTrapezoid('x', systmp, 'Area', nx*deltak(1) + maxBlipArea), commonRasterTime, systmp);
else
    gro = arg.gro;
end

% break up readout trapezoid
% Piece 1: ramp from 0 to end of gy/gz blip
% Piece 2: ADC window (symmetric)
% Piece 3: ramp from start of gy/gz blip to 0
[gro1, gro23] = mr.splitGradientAt(gro, blipDuration/2);
[gro2, gro3] = mr.splitGradientAt(gro23, mr.calcDuration(gro) - blipDuration/2);
gro2.delay = 0;
gro3.delay = 0;
gro1.delay = blipDuration/2;
gro31 = mr.addGradients({gro3, mr.scaleGrad(gro1, -1)}, sys);

% ADC event
% Number of readout samples must be multiple of 4 (TODO: check if this is actually needed)
if isempty(arg.adc) 
    TreadTmp = mr.calcDuration(gro2);
    numSamplesTmp = ceil(TreadTmp/dwell);
    numSamples = ceil(numSamplesTmp/4)*4;
    Tread = numSamples*dwell;
    adc = mr.makeAdc(numSamples, sys, ...
        'Duration', Tread, ...
        'Delay', blipDuration/2);  % + dwell/2? TODO
else
    adc = arg.adc;
end

% prephasers and spoilers
gxPre = trap4ge(mr.makeTrapezoid('x', sys, ...
    'Area', -gro.area/2), ...
    commonRasterTime, sys);
Tpre = mr.calcDuration(gxPre) + 4*commonRasterTime;
gyPre = trap4ge(mr.makeTrapezoid('y', sys, ...
    'Area', (kyInds(1)-ny/2)*deltak(2), ... 
    'Duration', Tpre-4*commonRasterTime), ...  % make a bit shorter than Tpre to ensure duration doesn't exceed Tpre after trap4ge
    commonRasterTime, sys);
if ~strcmp(tp, '3D')
    area = -floor(mb/2)*deltak(3);
else
    area = nz/2*deltak(3);
end
gzPre = trap4ge(mr.makeTrapezoid('z', sys, ...
    'Area', area, ...
    'Duration', Tpre-commonRasterTime), ...   % make < Tpre to ensure duration doesn't exceed Tpre
    commonRasterTime, sys);
gxSpoil = mr.makeTrapezoid('x', sys2, ...
    'Area', -nx*deltak(1)*nCyclesSpoil);
gzSpoil = mr.makeTrapezoid('z', sys2, ...
    'Area', nx*deltak(1)*nCyclesSpoil);


%% Calculate delays to achieve desired TE (for SMS and 3D) and TR.
% Note: gzRF includes gradient rephaser, while
% gzRF_se does NOT include gradient spoilers
kyIndAtTE = find(kyInds-ny/2 == min(abs(kyInds-ny/2)));
minTE = mr.calcDuration(gzRF) - mr.calcDuration(rf)/2 - rf.delay + mr.calcDuration(gxPre) + ...
        (kyIndAtTE-0.5) * mr.calcDuration(gro);
assert(TE+eps > minTE, sprintf('Requested TE < minimum TE (%f)', minTE));
TEdelay = floor((TE-minTE)/sys.blockDurationRaster) * sys.blockDurationRaster;

minTR = arg.fatSat*(mr.calcDuration(rfsat) + mr.calcDuration(gxSpoil)) + ...
    mr.calcDuration(gzRF) + TEdelay + ...
    mr.calcDuration(gxPre) + etl*mr.calcDuration(gro) + mr.calcDuration(gxSpoil);

if strcmp(tp, '3D')
    minTR = minTR + mr.calcDuration(gzPre);
end

TRdelay = round((TR/np-minTR-arg.segmentRingdownTime)/sys.blockDurationRaster) * sys.blockDurationRaster;
assert(TR > np*minTR, sprintf('Requested TR < minimum TR (%f)', minTR));

%% Slice/partition order
if strcmp(tp, '3D')
    IP = -nz/2:Rz:nz/2-1;
else
    % Interleaved slice ordering for SMS/2D
    % if even number of shots, swap last two to reduce slice cross-talk
    IP = [1:2:np 2:2:np];
    if ~mod(np,2)
        % for np = even, change order of last two partitions/shots
        l = length(IP);
        IP = IP([1:(l-2) l l-1]);
    end
end

%% Output trigger (stimulus trigger) event
trigOut = mr.makeDigitalOutputPulse('ext1', 'duration', 200e-6);

%% Assemble sequence
seq = mr.Sequence(sys);           

% temporal frame loop
rf_phase = 0;
rf_inc = 0;
msg = [];
for ifr = (1-nDummyFrames):nFrames

    for ii = 1:length(msg)
        fprintf('\b');
    end
    msg = sprintf('Frame %d of %d     ', ifr, nFrames);
    fprintf(msg);

    % dummy shots before turning on ADC, to reach steady state 
    isDummyShot = ifr < 1;

    % Segment/TR ID (see 'Pulseq on GE' manual)
    trid = 2 - isDummyShot;

    yBlipsOn = ~isDummyShot - eps; % trick: subtract eps to avoid scaling exactly to zero, while keeping scaling <1
    zBlipsOn = yBlipsOn - eps; 

    % slice (partition/SMS group) loop
    for p = IP
        rf.freqOffset = (1-strcmp(tp, '3D')) * round((p-1)*freq);  % frequency offset (Hz) for SMS slice shift

        % Label the start of segment instance
        seq.addBlock(mr.makeLabel('SET', 'TRID', trid));

        % fat sat and RF spoiling
        if arg.fatSat
            rfsat.phaseOffset = rf_phase/180*pi;
            seq.addBlock(rfsat);
            seq.addBlock(gxSpoil, gzSpoil);
            if arg.RFspoil
                rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
                rf_phase = mod(rf_phase+rf_inc, 360.0);
            end
        end

        % excitation pulse and RF spoiling
        rf.phaseOffset = rf_phase/180*pi - 2*pi*rf.freqOffset*mr.calcRfCenter(rf);  % align the phase for off-center slices
        adc.phaseOffset = rf_phase/180*pi;
        if ~arg.doNoiseScan
            seq.addBlock(rf, gzRF);
        else
            seq.addBlock(gzRF);
        end
        if arg.RFspoil
            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);
        end

        % TE delay
        if arg.trigOut
            seq.addBlock(mr.makeDelay(TEdelay), trigOut);
        else
            seq.addBlock(mr.makeDelay(TEdelay));
        end

        % Readout pre-phasers
        if strcmp(tp, '3D')
            seq.addBlock(gxPre, mr.scaleGrad(gyPre, arg.gySign*yBlipsOn), mr.scaleGrad(gzPre, p/(nz/2)*zBlipsOn*arg.gzPreOn));
        else
            seq.addBlock(gxPre, mr.scaleGrad(gyPre, arg.gySign*yBlipsOn), mr.scaleGrad(gzPre, zBlipsOn));
        end

        % first echo
        seq.addBlock(gro1);
        seq.addBlock(gro2, adc);

        % middle echoes
        for e = 2:(etl-1)
            seq.addBlock(mr.scaleGrad(gro23, (-1)^(e)), ...
                mr.scaleGrad(gyBlip, arg.gySign*yBlipsOn*kyStep(e-1)/max(kyStepMax,1)), ...
                mr.scaleGrad(gzBlip, zBlipsOn*kzStep(e-1)/max(kzStepMax,1)));
            seq.addBlock(mr.scaleGrad(gro2, (-1)^(e+1)), adc);
        end

        % last echo
        seq.addBlock(mr.scaleGrad(gro23, (-1)^(e+1)), ...
            mr.scaleGrad(gyBlip, arg.gySign*yBlipsOn*kyStep(e-1)/max(kyStepMax,1)), ...
            mr.scaleGrad(gzBlip, zBlipsOn*kzStep(e-1)/max(kzStepMax,1))); 
        seq.addBlock(mr.scaleGrad(gro2, (-1)^(e+2)), adc);
        seq.addBlock(mr.scaleGrad(gro1, (-1)^(e+2)));
        end

        % gz rephaser for 3D
        if strcmp(tp, '3D')
            seq.addBlock(mr.scaleGrad(gzPre, -p/(nz/2)*zBlipsOn*arg.gzPreOn));
        end

        % spoil. Disabling makes the k-space plot neater
        seq.addBlock(gxSpoil, gzSpoil);

        % TR delay
        seq.addBlock(mr.makeDelay(TRdelay));

        keyboard
    end
end
fprintf('\n');

%% If noise scan, add dummy rf pulse at the end since seg2ceq()
%% requires at least one rf pulse to be present in sequence
if arg.doNoiseScan
    seq.addBlock(mr.makeLabel('SET', 'TRID', 3*trid));
    seq.addBlock(rf, gzRF, mr.makeDelay(1));
end

%% Check sequence timing
[ok, error_report]=seq.checkTiming;
if (ok)
    fprintf('Timing check passed\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end


%% Write .seq file
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', arg.seqName);
ifn = [arg.seqName '.seq'];
seq.write(ifn);       % Write to pulseq file

% seq.plot('timeRange', [0 0.06]);

if ~arg.toGE
    return;
end

% add caipi.mat to the .tar file
%system(sprintf('tar --append --file=%s caipi.mat', ofn));

%% Optional slow step, but useful for testing during development,
%% e.g., for the real TE, TR or for staying within slewrate limits
% rep = seq.testReport;
% fprintf([rep{:}]);
