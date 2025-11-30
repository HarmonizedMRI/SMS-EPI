function [gro, adc] = writeEPI(voxelSize, N, TE, TR, alpha, mb, IY, IZ, nFrames, type, varargin)
% function [gro, adc] = writeEPI(voxelSize, N, TE, TR, alpha, mb, pf_ky, Ry, Rz, nFrames, type, varargin)
%
% SMS-EPI sequence in Pulseq
%
% Inputs:
%   voxelSize    [3]      meter
%   N            [3]      image matrix size 
%   TE           [1]      sec
%   alpha        [1]      flip angle (degrees)
%   mb           [1]      multiband/SMS factor
%   IY       [etl]    ky phase encoding indices, range is 1...ny, centered at ny/2+1
%   IZ       [etl]    kz partition encoding indices, range is 1...mb, centered at mb/2+1
%   Rz           [1]      kz undersampling factor
%   nFrames      [1]      number of time frames (image volumes)
%   type         string   'SMS' or '3D'
%
% Keyword-argument input options:
%   gro          struct   readout gradient struct. Default: create a new one and return from this function
%   adc          struct   acquisition struct. Default: create a new one and return from this function
%
% Outputs:
%   gro 
%   adc

% copy/modify inputs as needed
[nx ny nz] = deal(N(1), N(2), N(3));

is3D = strcmp(type, '3D');

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
arg.segmentRingdownTime = 117e-6;    % segment ringdown time for Pulseq on GE 
arg.doRefScan = false;               % do EPI ghost reference scan in frame 1
arg.doNoiseScan = false;
arg.seqName = 'epi';
arg.gySign = +1;
arg.freqSign = +1;
arg.fatFreqSign = -1;
arg.doConj = false;
arg.trigOut = false;

arg = toppe.utils.vararg_pair(arg, varargin);

arg.fatSat = pge2.utils.iff(arg.doNoiseScan, false, arg.fatSat);

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

slThick = pge2.utils.iff(is3D, 0.85*fov(3), fov(3)/nz);

dwell = 4e-6;                    % ADC sample time (s)

nCyclesSpoil = 2;    % number of spoiler cycles, along x and z
rfSpoilingInc = pge2.utils.iff(arg.RFspoil, 117, 0);    % RF spoiling increment (degrees)

rfTB = pge2.utils.iff(is3D, 8, 6);   % RF pulse time-bandwidth product
rfDur = pge2.utils.iff(is3D, 4e-3, 8e-3);  % RF pulse duration (s)

fatChemShift = 3.5*1e-6;                        % 3.5 ppm
fatOffresFreq = sys.gamma*sys.B0*fatChemShift;  % Hz


%% Create fat sat pulse 

fatsat.flip    = 90;      % degrees
fatsat.slThick = 1e5;     % dummy value (determines slice-select gradient, but we won't use it; just needs to be large to reduce dead time before+after rf pulse)
fatsat.tbw     = 3.5;     % time-bandwidth product
fatsat.dur     = 8.0;     % pulse duration (ms)

% RF waveform in Gauss
wav = toppe.utils.rf.makeslr(fatsat.flip, fatsat.slThick, fatsat.tbw, fatsat.dur, 1e-6, toppe.systemspecs(), ...
    'type', 'ex', ...    % fatsat pulse is a 90 so is of type 'ex', not 'st' (small-tip)
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


%% Create SMS excitation pulse
sysGE = toppe.systemspecs('maxGrad', sys.maxGrad/sys.gamma*100, ...   % G/cm
    'maxSlew', 0.7*sys.maxSlew/sys.gamma/10, ...           % G/cm/ms
    'maxRF', 0.15);
sliceSep = fov(3)/mb;   % center-to-center separation between SMS slices (m)
% freq = frequency offset (Hz) corresponding to a sliceSep
[rf, gzRF, freq, t_rf_center] = getsmspulse(alpha, slThick, rfTB, rfDur, ...
    mb, sliceSep, sysGE, sys, ...
    'doSim', arg.simulateSliceProfile, ...    % Plot simulated SMS slice profile
    'type', 'st', ...     % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
    'noRfOffset', is3D, ...   % don't shift slice (slab) for 3D
    'ftype', pge2.utils.iff(is3D, 'min', 'ls'));   % filter design. 'ls' = least squares, 'min' = minimum phase

rf.use = 'excitation';

freq = arg.freqSign * freq;

rf.signal = pge2.utils.iff(arg.doConj, conj(rf.signal), rf.signal);
rf.signal = pge2.utils.iff(arg.doNoiseScan, 1e-9 * rf.signal, rf.signal);  % practically zero


%% ky/kz encoding blip amplitude along echo train (multiples of deltak)
[IYlabel, kyStep] = ky2blipsandrewinders(IY);
kyStepMax = max(abs(kyStep(IYlabel)));
kyRewindMax = max(abs(kyStep(~IYlabel)));
kzStep = diff(IZ);
kzStepMax = max(abs(kzStep));

etl = length(IY);


%% Define readout gradients and ADC event

deltak = 1./fov;

commonRasterTime = 20e-6;   

% Start with the blips
gyBlip = mr.makeTrapezoid('y', sys, 'Area', kyStepMax*deltak(2));
gzBlip = mr.makeTrapezoid('z', sys, 'Area', kzStepMax*deltak(3));

blipDuration = max(mr.calcDuration(gyBlip), mr.calcDuration(gzBlip));
maxBlipArea = max(gyBlip.area, gzBlip.area);

% y rewinder
gyRewind = mr.makeTrapezoid('y', sys', 'Area', kyRewindMax*deltak(2));

% Readout trapezoid
if isempty(arg.gro) 
    %systmp = sys;
    %systmp.maxGrad = deltak(1)/dwell;  % to ensure >= Nyquist sampling
    %gro = mr.makeTrapezoid('x', systmp, 'Area', nx*deltak(1) + maxBlipArea); % this can fail for some reason
    gro = pge2.utils.makeTrapezoid('x', sys, 'Area', nx*deltak(1) + maxBlipArea, ...
        'maxGrad', deltak(1)/dwell);  
else
    gro = arg.gro;
end

% Break up readout trapezoid into parts.
% This way, we don't have to split the gy/gz blips,
% which in turn allows them to be arbitrarily scaled in the scan loop
% without having to store unique shapes for each scaling factor.
%
% gro_t0_t1: ramp from 0 to end of gy/gz blip (t1)
% gro_t1_t4: from t1 to t4. Contains ADC window (symmetric).
% gro_t4_t5: ramp from start of gy/gz blip (t4) to 0 
%
%   ^
%   |                         
%   |            +-------------+  gro.amplitude
%   |           /               \
%   |          /                 \
%   |         /                   \   
%   |  gamp  +                     +  gamp
%   |       /                       \
%   |      /                         \
%   +-----+---------------------------+--+-------> time
%         t0 t1 t2             t3  t4 t5 t6

gamp = gro.amplitude/gro.riseTime*blipDuration/2;

t0 = 0;
t1 = t0 + blipDuration/2;
t2 = t0 + gro.riseTime;
t3 = t2 + gro.flatTime;
t5 = t3 + gro.fallTime;
t4 = t5 - blipDuration/2;
t6 = t5 + blipDuration/2;

area_t0_t1 = gamp * (t1-t0)/2;
area_t1_t2 = (gro.amplitude-gamp) * (t2-t1)/2;
area_t2_t3 = (gro.amplitude-gamp) * (t3-t2);
area_t3_t4 = area_t1_t2;
area_t4_t5 = area_t0_t1;
area_t5_t6 = -area_t0_t1;

gro_t0_t1 = struct('type', 'grad', ...
      'channel', 'x', ...
      'waveform', [0 gamp], ...
      'delay', 0, ...
      'tt', [t0 t1] - t0, ...
      'shape_dur', t1 - t0, ...
      'area', area_t0_t1, ...
      'first', 0, ...
      'last', gamp);

gro_t1_t5 = struct('type', 'grad', ...
      'channel', 'x', ...
      'waveform', [gamp gro.amplitude gro.amplitude 0], ...
      'delay', 0, ...
      'tt', [t1 t2 t3 t5] - t1, ...
      'shape_dur', t5 - t1, ...
      'area', area_t1_t2 + area_t2_t3 + area_t3_t4 + area_t4_t5, ...
      'first', gamp, ...
      'last', 0);

gro_t1_t6 = struct('type', 'grad', ...
      'channel', 'x', ...
      'waveform', [gamp gro.amplitude gro.amplitude -gamp], ...
      'delay', 0, ...
      'tt', [t1 t2 t3 t6] - t1, ...
      'shape_dur', t6 - t1, ...
      'area', gro_t1_t5.area + area_t5_t6, ...
      'first', gamp, ...
      'last', -gamp);

% set blip delays (they happen at end of t1-t6 block)
gyBlip.delay = t5 - mr.calcDuration(gyBlip)/2 - t1;
gzBlip.delay = t5 - mr.calcDuration(gzBlip)/2 - t1;

% ADC event 
if isempty(arg.adc) 
    numSamples = sys.adcSamplesDivisor*round((t4-t1)/dwell/sys.adcSamplesDivisor);
    adc = mr.makeAdc(numSamples, sys, ...
        'Duration', dwell*numSamples, ...
        'Delay', 0);
else
    adc = arg.adc;
end

% prephasers and spoilers
gxPre = mr.makeTrapezoid('x', sys, ...
    'Area', -gro.area/2);
Tpre = mr.calcDuration(gxPre) + 4*commonRasterTime;
gyPre = mr.makeTrapezoid('y', sys, ...
    'Area', (IY(1)-ny/2-1)*deltak(2), ... 
    'Duration', Tpre-4*commonRasterTime);   % make a bit shorter than Tpre to ensure duration doesn't exceed Tpre after trap4ge
gzPre = mr.makeTrapezoid('z', sys, ...
    'Area', pge2.utils.iff(is3D, nz/2*deltak(3), (IZ(1)-mb/2-1)*deltak(3)), ...
    'Duration', Tpre-commonRasterTime);    % make < Tpre to ensure duration doesn't exceed Tpre
gxSpoil = mr.makeTrapezoid('x', sys2, ...
    'Area', -nx*deltak(1)*nCyclesSpoil);
gzSpoil = mr.makeTrapezoid('z', sys2, ...
    'Area', nx*deltak(1)*nCyclesSpoil);


%% Calculate delays to achieve desired TE and TR.
if 0
kyIndAtTE = find(IY-ny/2-1 == min(abs(IY-ny/2-1)));
minTE = mr.calcDuration(gzRF) - mr.calcDuration(rf)/2 - rf.delay + mr.calcDuration(gxPre) + ...
        (kyIndAtTE-0.5) * mr.calcDuration(gro);
assert(TE+eps > minTE, sprintf('Requested TE < minimum TE (%f)', minTE));
TEdelay = floor((TE-minTE)/sys.blockDurationRaster) * sys.blockDurationRaster;

minTR = arg.fatSat*(mr.calcDuration(rfsat) + mr.calcDuration(gxSpoil)) + ...
    mr.calcDuration(gzRF) + TEdelay + ...
    mr.calcDuration(gxPre) + etl*mr.calcDuration(gro) + mr.calcDuration(gxSpoil);

minTR = minTR + pge2.utils.iff(is3D, mr.calcDuration(gzPre), 0);

TRdelay = round((TR/np-minTR-arg.segmentRingdownTime)/sys.blockDurationRaster) * sys.blockDurationRaster;
assert(TR > np*minTR, sprintf('Requested TR < minimum TR (%f)', minTR));
end

%% Slice/partition order
if is3D
    IP = -nz/2:Rz:nz/2-1;
else
    % Interleaved slice ordering for SMS/2D
    IP = [1:2:np 2:2:np];
    if mod(np,2) == 0
        % for np = even, change order of last two partitions/shots to reduce slice cross-talk
        IP([end-1 end]) = IP([end end-1]);
    end
end

%% Output trigger (stimulus trigger) event
trigOut = pge2.utils.iff(arg.trigOut, mr.makeDigitalOutputPulse('ext1', 'duration', 200e-6), []);

%% Assemble sequence
seq = mr.Sequence(sys);           

% temporal frame loop
rf_phase = 0;
rf_inc = 0;

for ifr = 1:nFrames
    fprintf('\rFrame %d of %d     ', ifr, nFrames);

    yBlipsOn = ~arg.doRefScan - eps; % trick: subtract eps to avoid scaling exactly to zero, while keeping scaling <1
    zBlipsOn = yBlipsOn - eps; 

    % slice (partition/SMS group) loop
    for p = IP

        % Label the start of segment instance
        seq.addBlock(mr.makeLabel('SET', 'TRID', 1));

        % fat sat
        if arg.fatSat
            rfsat.phaseOffset = rf_phase/180*pi;
            seq.addBlock(rfsat);
            seq.addBlock(gxSpoil, gzSpoil);
            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);
        end

        % SMS excitation 
        rf.freqOffset = pge2.utils.iff(is3D, 0, round((p-1)*freq));  
        rf.phaseOffset = rf_phase/180*pi - 2*pi*rf.freqOffset*mr.calcRfCenter(rf);  % align the phase for off-center slices
        adc.phaseOffset = rf_phase/180*pi;
        seq.addBlock(rf, gzRF);
        rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase = mod(rf_phase+rf_inc, 360.0);

        % TE delay and trigger output pulse
        %seq.addBlock(mr.makeDelay(TEdelay), trigOut);

        % Readout pre-phasers
        amp = pge2.utils.iff(is3D, p/(nz/2)*zBlipsOn*arg.gzPreOn, zBlipsOn);
        seq.addBlock(gxPre, mr.scaleGrad(gyPre, arg.gySign*yBlipsOn), mr.scaleGrad(gzPre, amp));

        % echo train
        seq.addBlock(gro_t0_t1);
        
        for e = 1:etl-1
            if ~IYlabel(e)           % ky rewinder
                seq.addBlock(mr.scaleGrad(gro_t1_t5, (-1)^(e+1)), adc);
                tmpDelay = gzBlip.delay;
                gzBlip.delay = 0;
                seq.addBlock(mr.scaleGrad(gyRewind, arg.gySign*yBlipsOn*kyStep(e)/max(kyRewindMax,1)), ...
                    mr.scaleGrad(gzBlip, zBlipsOn*kzStep(e)/max(kzStepMax,1)));
                gzBlip.delay = tmpDelay;
                seq.addBlock(mr.scaleGrad(gro_t0_t1, (-1)^(e+2)));
            else                     % ky blip
                seq.addBlock(adc, ...
                    mr.scaleGrad(gro_t1_t6, (-1)^(e+1)), ...
                    mr.scaleGrad(gyBlip, arg.gySign*yBlipsOn*kyStep(e)/max(kyStepMax,1)), ...
                    mr.scaleGrad(gzBlip, zBlipsOn*kzStep(e)/max(kzStepMax,1)));
            end
        end

        seq.addBlock(mr.scaleGrad(gro_t1_t5, (-1)^(etl+1)), adc);

        % finish out the TR
        if is3D
            seq.addBlock(mr.scaleGrad(gzPre, -amp));
        end
        seq.addBlock(gxSpoil, gzSpoil);
        %seq.addBlock(mr.makeDelay(TRdelay));
    end
end
fprintf('\n');

% Noise scan (add gaps to make room for adc dead/ringdown time)
if arg.doNoiseScan
    seq.addBlock(mr.makeLabel('SET', 'TRID', 2));
    seq.addBlock(mr.makeDelay(1));
    for ii = 1:10
        seq.addBlock(adc);
        seq.addBlock(mr.makeDelay(0.2));
    end
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

seq.plot('timeRange', [0 0.1], 'stacked', true);

% add caipi.mat to the .tar file
%system(sprintf('tar --append --file=%s caipi.mat', ofn));

%% Optional slow step, but useful for testing during development,
%% e.g., for the real TE, TR or for staying within slewrate limits
% rep = seq.testReport;
% fprintf([rep{:}]);
