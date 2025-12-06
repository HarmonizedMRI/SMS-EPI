function seq = writeEPI(seqName, sys, voxelSize, N, TR, alpha, mb, IY, IZ, nFrames, type, opts)
% Build an SMS/3D EPI Pulseq sequence.
%   seq = writeEPI(seqName, sys, voxelSize, N, TR, alpha, mb, IY, IZ, ...
%                  nFrames, type, opts)
%
% Inputs
%   seqName     string or []    output file name (with or without .seq extension)
%   sys         struct          scanner info, see mr.opts()
%   voxelSize   [3]             voxel dimensions [m]
%   N           [3]             image matrix size 
%   TR          [1] or 'min'    specify volume TR [s] or use minimum
%   alpha       [1]             flip angle (degrees)
%   mb          [1]             multiband/SMS factor
%   IY          [etl]           ky phase encoding indices, range is 1...ny centered at ny/2+1
%   IZ          [etl]           kz partition encoding indices, range is 1...mb centered at mb/2+1
%   nFrames     [1]             number of time frames (image volumes)
%   type        string          'SMS' or '3D'
%   opts        optional struct with fields that override defaults

% Input checks
narginchk(11,12);
if nargin < 12 || isempty(opts), opts = struct(); end

validateattributes(voxelSize, {'numeric'}, {'vector','numel',3});
validateattributes(N, {'numeric'}, {'vector','integer','numel',3});
validateattributes(alpha, {'numeric'}, {'scalar'});
validateattributes(mb, {'numeric'}, {'scalar','positive','integer'});
validateattributes(IY, {'numeric'}, {'vector','integer'});
validateattributes(IZ, {'numeric'}, {'vector','integer'});
validateattributes(nFrames, {'numeric'}, {'scalar','integer','positive'});
assert(ismember(type, {'SMS','3D'}), 'type must be ''SMS'' or ''3D''.');

assert(length(IY) == length(IZ), 'IY and IZ must be the same length');

% local namespace struct
lv = struct();

[lv.nx lv.ny lv.nz] = deal(N(1), N(2), N(3));

lv.is3D = strcmp(type, '3D');

assert(mod(lv.nz, mb) == 0, 'nz must be divisible by mb.');
lv.np = lv.nz / mb;   % number of partitions (excitations per volume)

fprintf('mb=%d\n', mb); 

% --- Defaults ---
arg.fatSat = true;  
arg.RFspoil = true;
arg.gySign = +1;
arg.freqSign = +1;
arg.fatFreqSign = -1;
arg.doConj = false;
arg.doRefScan = false;          % turn off ky and ky encoding
arg.trigOut = false;
arg.doNoiseScan = false;
arg.TEdelay = 0;
arg.segmentRingdownTime = 117e-6;  % segment ringdown time for Pulseq on GE 
arg.simulateSliceProfile = false;  % simulate SMS profile and display
arg.gzPreOn = true;                % prephase along kz by -mb/2*deltak
arg.plot = false;

% Overwrite defaults with user-specified fields
fn = fieldnames(opts);
for k = 1:numel(fn)
    arg.(fn{k}) = opts.(fn{k});
end

% --- Some more parameters ---
fov = voxelSize .* [lv.nx lv.ny lv.nz];       % FOV (m)

slThick = pge2.utils.iff(lv.is3D, 0.85*fov(3), fov(3)/lv.nz);

dwell = 4e-6;                    % ADC sample time (s)

nCyclesSpoil = 2;    % number of spoiler cycles, along x and z
lv.rfSpoilingInc = pge2.utils.iff(arg.RFspoil, 117, 0);    % RF spoiling increment (degrees)

rfTB = pge2.utils.iff(lv.is3D, 8, 6);   % RF pulse time-bandwidth product
rfDur = pge2.utils.iff(lv.is3D, 4e-3, 8e-3);  % RF pulse duration (s)

fatChemShift = 3.5*1e-6;                        % 3.5 ppm
fatOffresFreq = sys.gamma*sys.B0*fatChemShift;  % Hz


% --- Create fat sat pulse ---

fatsat.flip    = 90;      % degrees
fatsat.slThick = 1e5;     % dummy value (determines slice-select gradient, but we won't use it; just needs to be large to reduce dead time before+after rf pulse)
fatsat.tbw     = 3.5;     % time-bandwidth product
fatsat.dur     = 8.0;     % pulse duration (ms)

% RF waveform in Gauss. Designed on 4us raster.
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
lv.rfsat = mr.makeArbitraryRf(rfp, ...
    flip*abs(sum(rfp*sys.rfRasterTime))*(2*pi), ...
    'use', 'excitation', ...
    'system', sys);
lv.rfsat.signal = lv.rfsat.signal/max(abs(lv.rfsat.signal))*max(abs(rfp)); % ensure correct amplitude (Hz)
lv.rfsat.freqOffset = arg.fatFreqSign*425;  % Hz

% --- SMS excitation pulse ---
sysGE = toppe.systemspecs('maxGrad', sys.maxGrad/sys.gamma*100, ...   % G/cm
    'maxSlew', 0.7*sys.maxSlew/sys.gamma/10, ...           % G/cm/ms
    'maxRF', 0.15);
sliceSep = fov(3)/mb;   % center-to-center separation between SMS slices (m)
% freq = frequency offset (Hz) corresponding to a sliceSep
[lv.rf, lv.gzRF, freq, t_rf_center] = hmriutils.rf.getsmspulse(alpha, slThick, rfTB, rfDur, ...
    mb, sliceSep, sysGE, sys, ...
    'doSim', arg.simulateSliceProfile, ...     % Plot simulated SMS slice profile
    'type', 'st', ...                          % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
    'noRfOffset', lv.is3D, ...                 % don't shift slice (slab) for 3D
    'ftype', pge2.utils.iff(lv.is3D, ...
    'min', 'ls'));                          % 'ls' = least squares, 'min' = minimum phase

% rf design is on 4us raster time (both gradient and RF) so now we interpolate
tt = sys.rfRasterTime/2:sys.rfRasterTime:lv.rf.shape_dur;
lv.rf.signal = interp1(lv.rf.t, lv.rf.signal, tt, 'linear', 'extrap');
lv.rf.tt = tt;

% make sure the gradient delay is on raster boundary
dpad = sys.gradRasterTime - rem(lv.gzRF.delay, sys.gradRasterTime);
lv.gzRF.delay = lv.gzRF.delay + dpad;
lv.rf.delay = lv.rf.delay + dpad;

% some final settings
lv.rf.use = 'excitation';
lv.freq = arg.freqSign * freq;
lv.rf.signal = pge2.utils.iff(arg.doConj, conj(lv.rf.signal), lv.rf.signal);

% --- ky/kz encoding blip amplitude along echo train (multiples of deltak) ---
[lv.IYlabel, lv.kyStep] = ky2blipsandrewinders(IY);
lv.kyStepMax = max(abs(lv.kyStep(lv.IYlabel)));
if sum(lv.IYlabel) == length(IY) - 1   % no rewinders are present
    lv.kyRewindMax = 1e-9;  % avoid exactly zero
else
    lv.kyRewindMax = max(abs(lv.kyStep(~lv.IYlabel)));
end    
lv.kzStep = diff(IZ);
lv.kzStepMax = max(abs(lv.kzStep)) + 1e-9;  % avoid exactly zero

lv.etl = length(IY);

% --- Readout gradients and ADC event ---

deltak = 1./fov;

% Start with the blips
lv.gyBlip = mr.makeTrapezoid('y', sys, 'Area', lv.kyStepMax*deltak(2));
lv.gzBlip = mr.makeTrapezoid('z', sys, 'Area', lv.kzStepMax*deltak(3));

blipDuration = max(mr.calcDuration(lv.gyBlip), mr.calcDuration(lv.gzBlip));
maxBlipArea = max(lv.gyBlip.area, lv.gzBlip.area);

% y rewinder
lv.gyRewind = mr.makeTrapezoid('y', ...
    pge2.utils.setfields(sys, 'maxSlew', sys.maxSlew*0.8), ...
    'Area', lv.kyRewindMax*deltak(2));

% Readout trapezoid
lv.gro = pge2.utils.makeTrapezoid('x', sys, 'Area', lv.nx*deltak(1) + maxBlipArea, ...
    'maxGrad', deltak(1)/dwell);  

% Split readout into shape parts to avoid unique shapes per block.
%
% gro_t0_t1: ramp from 0 to end of gy/gz blip (t1)
% gro_t1_t4: from t1 to t4. Contains ADC window (symmetric).
% gro_t4_t5: ramp from start of gy/gz blip (t4) to 0 
%
%   ^
%   |            +-------------+   gro.amplitude
%   |           /               \
%   |          /                 \
%   |         /                   \   
%   |        +<------- ADC ------->+   gamp
%   |       /                       \
%   |      /                         \
%   +-----+---------------------------+--+-------> time
%         t0 t1 t2             t3  t4 t5 t6

gamp = lv.gro.amplitude/lv.gro.riseTime*blipDuration/2;

t0 = 0;
t1 = t0 + blipDuration/2;
t2 = t0 + lv.gro.riseTime;
t3 = t2 + lv.gro.flatTime;
t5 = t3 + lv.gro.fallTime;
t4 = t5 - blipDuration/2;
t6 = t5 + blipDuration/2;

area_t0_t1 = gamp * (t1-t0)/2;
area_t1_t2 = (lv.gro.amplitude-gamp) * (t2-t1)/2;
area_t2_t3 = (lv.gro.amplitude-gamp) * (t3-t2);
area_t3_t4 = area_t1_t2;
area_t4_t5 = area_t0_t1;
area_t5_t6 = -area_t0_t1;

lv.gro_t0_t1 = struct('type', 'grad', ...
      'channel', 'x', ...
      'waveform', [0 gamp], ...
      'delay', 0, ...
      'tt', [t0 t1] - t0, ...
      'shape_dur', t1 - t0, ...
      'area', area_t0_t1, ...
      'first', 0, ...
      'last', gamp);

lv.gro_t1_t5 = struct('type', 'grad', ...
      'channel', 'x', ...
      'waveform', [gamp lv.gro.amplitude lv.gro.amplitude 0], ...
      'delay', 0, ...
      'tt', [t1 t2 t3 t5] - t1, ...
      'shape_dur', t5 - t1, ...
      'area', area_t1_t2 + area_t2_t3 + area_t3_t4 + area_t4_t5, ...
      'first', gamp, ...
      'last', 0);

lv.gro_t1_t6 = struct('type', 'grad', ...
      'channel', 'x', ...
      'waveform', [gamp lv.gro.amplitude lv.gro.amplitude -gamp], ...
      'delay', 0, ...
      'tt', [t1 t2 t3 t6] - t1, ...
      'shape_dur', t6 - t1, ...
      'area', lv.gro_t1_t5.area + area_t5_t6, ...
      'first', gamp, ...
      'last', -gamp);

% set blip delays (they happen at end of t1-t6 block)
lv.gyBlip.delay = t5 - mr.calcDuration(lv.gyBlip)/2 - t1;
lv.gzBlip.delay = t5 - mr.calcDuration(lv.gzBlip)/2 - t1;

% ADC event 
numSamples = sys.adcSamplesDivisor*round((t4-t1)/dwell/sys.adcSamplesDivisor);
lv.adc = mr.makeAdc(numSamples, sys, ...
    'Duration', dwell*numSamples, ...
    'Delay', 0);

% prephasers. Reduce slew a bit.
lv.gxPre = mr.makeTrapezoid('x', ...
    pge2.utils.setfields(sys, 'maxSlew', sys.maxSlew*0.8), ...
    'Area', -lv.gro.area/2);
Tpre = mr.calcDuration(lv.gxPre);
lv.gyPre = mr.makeTrapezoid('y', ...
    pge2.utils.setfields(sys, 'maxSlew', sys.maxSlew*0.8), ...
    'Area', (IY(1)-lv.ny/2-1)*deltak(2), ... 
    'Duration', Tpre);
lv.gzPre = mr.makeTrapezoid('z', ...
    pge2.utils.setfields(sys, 'maxSlew', sys.maxSlew*0.8), ...
    'Area', pge2.utils.iff(lv.is3D, lv.nz/2*deltak(3), (IZ(1)-mb/2-1)*deltak(3)), ...
    'Duration', Tpre); 

% spoilers. Reduce slew a lot.
lv.gxSpoil = mr.makeTrapezoid('x', ...
    pge2.utils.setfields(sys, 'maxSlew', sys.maxSlew*0.5), ...
    'Area', -lv.nx*deltak(1)*nCyclesSpoil);
lv.gzSpoil = mr.makeTrapezoid('z', ...
    pge2.utils.setfields(sys, 'maxSlew', sys.maxSlew*0.5), ...
    'Area', lv.nx*deltak(1)*nCyclesSpoil);

% --- Slice/partition order ---
if lv.is3D
    IP = -lv.nz/2:Rz:lv.nz/2-1;
else
    % Interleaved slice ordering for SMS/2D
    IP = [1:2:lv.np 2:2:lv.np];
    if mod(lv.np,2) == 0
        % for lv.np = even, change order of last two partitions/shots to reduce slice cross-talk
        IP([end-1 end]) = IP([end end-1]);
    end
end

% --- Output trigger (stimulus trigger) event ---
trigOut = pge2.utils.iff(arg.trigOut, mr.makeDigitalOutputPulse('ext1', 'duration', 200e-6), []);

% --- Assemble sequence ---
seq = mr.Sequence(sys);           

rf_phase = 0;
rf_inc = 0;

lv.yBlipsOn = ~arg.doRefScan - eps; % trick: subtract eps to avoid scaling exactly to zero, while keeping scaling <1
lv.zBlipsOn = lv.yBlipsOn; 

for ifr = 1:nFrames
    fprintf('\rFrame %d of %d     ', ifr, nFrames);

    % slice (partition/SMS group) loop
    for p = IP
        % Label the start of segment instance
        seq.addBlock(mr.makeLabel('SET', 'TRID', 1));

        % add a TR: fat sat => SMS slice excitation => EPI readout)
        [seq, rf_phase, rf_inc] = sub_addEPIshot(seq, lv, arg, p, rf_phase, rf_inc);

        % add delay to achieve requested TR
        if p == IP(1) 
            minTR = seq.duration + arg.segmentRingdownTime;
            if ischar(TR)
                TR = lv.np * minTR;
            end
            assert(TR >= lv.np*minTR, sprintf('Requested TR (%.3f ms) < minimum TR (%.3f ms)', TR, lv.np*minTR));
            TRdelay = round((TR/lv.np - minTR)/sys.blockDurationRaster) * sys.blockDurationRaster;
        end
        seq.addBlock(mr.makeDelay(TRdelay));
    end
end
fprintf('\n');

% If noise scan, add dummy RF pulse at end so the sequence contains at least one RF pulse
if arg.doNoiseScan
    seq.addBlock(mr.makeLabel('SET', 'TRID', 2));
    seq.addBlock(lv.rf, lv.gzRF, mr.makeDelay(0.1));
end

% --- Check sequence timing ---
[ok, error_report]=seq.checkTiming;
if (ok)
    fprintf('Timing check passed\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

% --- Write .seq file ---
seq.setDefinition('FOV', fov);
if ~isempty(seqName)
    seqName = replace(seqName, {'.seq', '.pge'}, '');
    seq.setDefinition('Name', seqName);
    seq.write([seqName '.seq']);
end

if arg.plot
    seq.plot('timeRange', [0 TR/lv.np], 'stacked', true);
end

% add caipi.mat to the .tar file
%system(sprintf('tar --append --file=%s caipi.mat', ofn));

% Optional slow step, but useful for testing during development,
% e.g., for the real TE, TR or for staying within slewrate limits
% rep = seq.testReport;
% fprintf([rep{:}]);



%% Subfunction: add one EPI shot
function [sq, rf_phase, rf_inc] = sub_addEPIshot(sq, lv, arg, p, rf_phase, rf_inc)

    % fat sat
    if arg.fatSat & ~arg.doNoiseScan
        lv.rfsat.phaseOffset = rf_phase/180*pi;
        sq.addBlock(lv.rfsat);
        sq.addBlock(lv.gxSpoil, lv.gzSpoil);
        rf_inc = mod(rf_inc+lv.rfSpoilingInc, 360.0);
        rf_phase = mod(rf_phase+rf_inc, 360.0);
    end

    % SMS excitation 
    lv.rf.freqOffset = pge2.utils.iff(lv.is3D, 0, round((p-1)*lv.freq));  
    lv.rf.phaseOffset = rf_phase/180*pi - 2*pi*lv.rf.freqOffset*mr.calcRfCenter(lv.rf);  % align the phase for off-center slices
    lv.adc.phaseOffset = rf_phase/180*pi;
    sq.addBlock(pge2.utils.iff(arg.doNoiseScan, [], lv.rf), lv.gzRF);

    rf_inc = mod(rf_inc+lv.rfSpoilingInc, 360.0);
    rf_phase = mod(rf_phase+rf_inc, 360.0);

    % TE delay and trigger output pulse
    %sq.addBlock(mr.makeDelay(TEdelay), trigOut);
    sq.addBlock(mr.makeDelay(arg.TEdelay));

    % Readout pre-phasers
    amp = pge2.utils.iff(lv.is3D, p/(lv.nz/2)*lv.zBlipsOn*arg.gzPreOn, lv.zBlipsOn);
    sq.addBlock(lv.gxPre, mr.scaleGrad(lv.gyPre, arg.gySign*lv.yBlipsOn), mr.scaleGrad(lv.gzPre, amp));

    % echo train
    sq.addBlock(lv.gro_t0_t1);

    for e = 1:lv.etl-1
        if lv.IYlabel(e)   % readout followed by ky/kz blip
            sq.addBlock(lv.adc, ...
                mr.scaleGrad(lv.gro_t1_t6, (-1)^(e+1)), ...
                mr.scaleGrad(lv.gyBlip, arg.gySign*lv.yBlipsOn*lv.kyStep(e)/max(lv.kyStepMax,1)), ...
                mr.scaleGrad(lv.gzBlip, lv.zBlipsOn*lv.kzStep(e)/max(lv.kzStepMax,1)));
        else               % readout followed by ky rewinder
            sq.addBlock(mr.scaleGrad(lv.gro_t1_t5, (-1)^(e+1)), lv.adc);
            tmpDelay = lv.gzBlip.delay;
            lv.gzBlip.delay = 0;
            sq.addBlock(mr.scaleGrad(lv.gyRewind, arg.gySign*lv.yBlipsOn*lv.kyStep(e)/max(lv.kyRewindMax,1)), ...
                mr.scaleGrad(lv.gzBlip, lv.zBlipsOn*lv.kzStep(e)/max(lv.kzStepMax,1)));
            lv.gzBlip.delay = tmpDelay;
            sq.addBlock(mr.scaleGrad(lv.gro_t0_t1, (-1)^(e+2)));
        end
    end

    sq.addBlock(mr.scaleGrad(lv.gro_t1_t5, (-1)^(lv.etl+1)), lv.adc);

    % finish out the TR
    if lv.is3D
        sq.addBlock(mr.scaleGrad(lv.gzPre, -amp));
    end

    sq.addBlock(lv.gxSpoil, lv.gzSpoil);

    return
