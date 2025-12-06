function writeB0(seqName, sys, voxelSize, N, alpha)
% function writeB0(seqName, sys, voxelSize, N, alpha)
%
% 3D GRE B0 mapping sequence 
%
% Inputs
%   seqName     string   .seq file name
%   sys         struct   see mr.opts()
%   voxelSize   [3]      m
%   N           [3]      matrix size
%   alpha       [1]      flip angle [degrees]

% --- Acquisition parameters ---
[nx ny nz] = deal(N(1), N(2), N(3));
assert(nz > 1, 'Only 3D scan is supported');
fov = voxelSize.*[nx ny nz];
dwell = 4*sys.adcRasterTime;    % ADC sample time (s)
fatChemShift = 3.5e-6;          % 3.5 ppm
fatOffresFreq = sys.gamma*sys.B0*fatChemShift;  % Hz
TE = 1/fatOffresFreq*[1 2];     % fat and water in phase for both echoes

rfSpoilingInc = 117;            % RF spoiling increment
nCyclesSpoil = 2;               % number of spoiler cycles
Tpre = 1.0e-3;                  % prephasing trapezoid duration

% --- Sequence elements ---

% non-selective pulse
[rf] = mr.makeBlockPulse(alpha/180*pi, sys, 'Duration', 0.2e-3, 'use', 'excitation');

% Define other gradients and ADC events
deltak = 1./fov;
Tread = nx*dwell;

gyPre = mr.makeTrapezoid('y', sys, ...
    'Area', ny*deltak(2)/2);   % PE1 gradient, max positive amplitude
gzPre = mr.makeTrapezoid('z', sys, ...
    'Area', nz*deltak(3)/2);   % PE2 gradient, max amplitude

gx = mr.makeTrapezoid('x', sys, ...  % readout trapezoid, temporary object
    'Amplitude', nx*deltak(1)/Tread, ...
    'FlatTime', Tread);
gxPre = mr.makeTrapezoid('x', sys, ...
    'Area', -gx.area/2);

adc = mr.makeAdc(nx, sys, ...
    'Duration', Tread,...
    'Delay', gx.riseTime);

gxSpoil = mr.makeTrapezoid('x', sys, ...
    'Area', nx*deltak(1)*nCyclesSpoil);
gzSpoil = mr.makeTrapezoid('z', sys, ...
    'Area', nx*deltak(1)*nCyclesSpoil);

% --- y/z PE steps. Avoid exactly zero ---
pe1Steps = ((0:ny-1)-ny/2)/ny*2 + 1e-9;
pe2Steps = ((0:nz-1)-nz/2)/nz*2 + 1e-9;

% --- Calculate timing ---
TEmin = rf.shape_dur/2 + rf.ringdownTime + mr.calcDuration(gxPre) ...
      + adc.delay + nx/2*dwell;
TEdelay = ceil((TE-TEmin)/sys.gradRasterTime)*sys.gradRasterTime;
if TEdelay < 0
    TE = 1/fatOffresFreq*[2 3]; 
    TEdelay = ceil((TE-TEmin)/sys.gradRasterTime)*sys.gradRasterTime;
end

% --- Loop over phase encodes and define sequence blocks ---

% iz < 0: Dummy shots to reach steady state (ADC off)
% iz = 0: ADC is turned on and used for receive gain calibration on GE scanners
% iz > 0: Image acquisition

seq = mr.Sequence(sys);           

nDummyZLoops = 1;
rf_phase = 0;
rf_inc = 0;

TRmin = 0;
%for iz = -nDummyZLoops:nz
for iz = 1:4
    isDummyTR = iz < 0;

    fprintf('\rz encode %d of %d ', iz, nz);

    for iY = 1:ny
        % Turn on/off y and z prephasing lobes
        % Avoid scaling to exactly zero so the shapes are preserved in the .seq file
        yStep = (iz > 0) * pe1Steps(iY) + eps;
        zStep = (iz > 0) * pe2Steps(max(1,iz)) + eps;

        for c = 1:length(TE)

            % update RF/adc phase (RF spoiling)
            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);
            rf.phaseOffset = rf_phase/180*pi;
            adc.phaseOffset = rf_phase/180*pi;

            % mark start of segment (block group) by adding TRID label
            seq.addBlock(mr.makeLabel('SET', 'TRID', 1 + isDummyTR));

            % excitation
            seq.addBlock(rf);
            
            % encoding
            seq.addBlock(mr.makeDelay(TEdelay(c)));   % variable delay block
            seq.addBlock(gxPre, ...
                mr.scaleGrad(gyPre, yStep), ...
                mr.scaleGrad(gzPre, zStep));
            if isDummyTR
                seq.addBlock(gx);
            else
                seq.addBlock(gx, adc);
            end

            % rephasing and spoiling
            seq.addBlock(gxSpoil, ...
                mr.scaleGrad(gyPre, -yStep), ...
                mr.scaleGrad(gzPre, -zStep));
            seq.addBlock(gzSpoil);

            % keep TR constant
            minVariableDelay = 8e-6;  % mimimum duration of variable delay block (WAIT pulse in EPIC)
            seq.addBlock(mr.makeDelay(minVariableDelay + TEdelay(end) - TEdelay(c)));

            TRmin = pge2.utils.iff(TRmin, TRmin, seq.duration); % for plotting below
        end
    end
end
fprintf('\n');

% --- Output for execution ---
pge2.utils.writeseq(seq, fov, seqName);

% --- Plot sequence ---
Noffset = length(TE)*ny*(nDummyZLoops+1);
seq.plot('timerange',[Noffset Noffset+2*length(TE)]*TRmin, 'timedisp', 'ms');

