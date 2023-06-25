% SMS-EPI sequence in Pulseq, optimized
% for identical execution on Siemens and GE scanners.
%
% This script creates the file 'smsepi.seq', that can be executed directly
% on Siemens MRI scanners using the Pulseq interpreter.
% The .seq file can also be converted to a .tar file that can be executed on GE
% scanners, see main.m.
%
% The experimental parameters below are chosen such that the sequence 
% can be executed identically (to us precision) on Siemens and GE systems.
% For more information about preparing a Pulseq file for execution on GE scanners,
% see https://github.com/jfnielsen/TOPPEpsdSourceCode/wiki.

%% Define experimental parameters
sys = mr.opts('maxGrad', 22, 'gradUnit','mT/m', ...
              'maxSlew', 80, 'slewUnit', 'T/m/s', ...
              'rfDeadTime', 100e-6, ...
              'rfRingdownTime', 60e-6, ...
              'adcDeadTime', 40e-6, ...
              'adcRasterTime', 2e-6, ...
              'gradRasterTime', 10e-6, ...
              'blockDurationRaster', 10e-6, ...
              'B0', 3.0);

timessi = 100e-6;    % start sequence interrupt (SSI) time (required delay at end of block group/TR)

voxelSize = [2.4 2.4 2.4]*1e-3;     % m
mb = 6;                             % multiband/SMS factor
Nx = 92; Ny = Nx; Nz = mb*10;       % Matrix size
fov = voxelSize .* [Nx Ny Nz];      % FOV (m)

dwell = 2e-6;                       % ADC sample time (s). For GE, must be multiple of 2us.
alpha = 60;                         % flip angle (degrees)

fatChemShift = 3.5*1e-6;            % 3.5 ppm
fatOffresFreq = sys.gamma*sys.B0*fatChemShift;  % Hz
TE = 30e-3;
nCyclesSpoil = 2;    % number of spoiler cycles, along x and z

rfTB  = 6;          % time-bandwidth product
rfDur = 8e-3;       % RF pulse duration (s)

%% Excitation pulse

% design pulse
clear rf ex
sysGE = toppe.systemspecs('maxGrad', sys.maxGrad/sys.gamma*100, ...   % G/cm
    'maxSlew', sys.maxSlew/sys.gamma/10, ...           % G/cm/ms
    'maxRF', 0.25);
sliceSep = fov(3)/mb;   % center-to-center separation between SMS slices (m)
[wav, g, freq] = getsmspulse(alpha, voxelSize(3), rfTB, rfDur, ...
    mb, sliceSep, sysGE, ...
    'doSim', false, ...   % Plot simulated SMS slice profile
    'type', 'st', ...     % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
    'ftype', 'ls');       % filter design. 'ls' = least squares

% Convert from Gauss (Gauss/cm) to Hz (Hz/m), and interpolate to sys.rf/gradRasterTime
% Pad with zeros to make equal length
wav = pulsegeq.rf2pulseq(wav, sysGE.raster, sys);  
g = pulsegeq.g2pulseq(g, sysGE.raster, sys);  
wavdur = length(wav)*sys.rfRasterTime;
gdur = length(g)*sys.gradRasterTime;
if gdur < wavdur
    n = ceil((wavdur-gdur)/sys.gradRasterTime);
    g = [g zeros(1, n)];
    gdur = length(g)*sys.gradRasterTime;
end
wav = [wav zeros(1, round((gdur-wavdur)/sys.rfRasterTime))];
rf = mr.makeArbitraryRf(wav, alpha/180*pi, ...
            'system', sys);

gzRf = mr.makeArbitraryGrad('z', g, sys);


return
% Convert from Gauss to Hz, and interpolate to sys.rfRasterTime


% Create a new sequence object
seq = mr.Sequence(sys);           

% Create non-selective pulse
[rf, rfDelay] = mr.makeBlockPulse(alpha/180*pi, sys, 'Duration', alphaPulseDuration);

% Define other gradients and ADC events
% Cut the redaout gradient into two parts for optimal spoiler timing
deltak = 1./fov;
Tread = Nx*dwell;
gxPre = mr.makeTrapezoid('x', sys, ...
    'Area', -Nx*deltak(1)/2, ...
    'Duration', Tpre);
gyPre = mr.makeTrapezoid('y', sys, ...
    'Area', Ny*deltak(2)/2, ...   % maximum PE1 gradient, max positive amplitude
    'Duration', Tpre);
gzPre = mr.makeTrapezoid('z', sys, ...
    'Area', Nz*deltak(3)/2, ...   % maximum PE2 gradient, max positive amplitude
    'Duration', Tpre);
gxtmp = mr.makeTrapezoid('x', sys, ...   % temporary object
    'Amplitude', Nx*deltak(1)/Tread, ...
    'FlatTime', Tread);
adc = mr.makeAdc(Nx, sys, ...
    'Duration', Tread,...
    'Delay', gxtmp.riseTime);
gxtmp2 = mr.makeTrapezoid('x', sys, ...  % temporary object
    'Amplitude', Nx*deltak(1)/Tread, ...
    'FlatTime', Tread + adc.deadTime);   % extend flat time so we can split at end of ADC dead time
gzSpoil = mr.makeTrapezoid('z', sys, ...
    'Area', Nx*deltak(1)*nCyclesSpoil);
gxSpoil = mr.makeExtendedTrapezoidArea('x', gxtmp.amplitude, 0, gzSpoil.area, sys);
[gx, ~] = mr.splitGradientAt(gxtmp2, gxtmp2.riseTime + gxtmp2.flatTime);

% y/z PE steps
pe1Steps = ((0:Ny-1)-Ny/2)/Ny*2;
pe2Steps = ((0:Nz-1)-Nz/2)/Nz*2;

% Calculate timing
TEmin = rf.shape_dur/2 + rf.ringdownTime + mr.calcDuration(gxPre) ...
      + adc.delay + Nx/2*dwell;
delayTE = ceil((TE-TEmin)/seq.gradRasterTime)*seq.gradRasterTime;
TRmin = mr.calcDuration(rf) + delayTE + mr.calcDuration(gxPre) ...
      + mr.calcDuration(gx) + mr.calcDuration(gxSpoil);
delayTR = ceil((TR-TRmin)/seq.gradRasterTime)*seq.gradRasterTime;

% Loop over phase encodes and define sequence blocks
% iZ < 0: Dummy shots to reach steady state
% iZ = 0: ADC is turned on and used for receive gain calibration on GE scanners (during auto prescan)
% iZ > 0: Image acquisition
nDummyZLoops = 2;
for iZ = -nDummyZLoops:Nz
    if iZ > 0
        for ib = 1:40
            fprintf('\b');
        end
        fprintf('Writing kz encode %d of %d', iZ, Nz);
    end
    for iY = 1:Ny
        % turn off y and z prephasing lobes during receive gain calibration (auto prescan)
        yStep = (iZ > 0) * pe1Steps(iY);
        zStep = (iZ > 0) * pe2Steps(max(1,iZ));
        for c = 1:length(TE)
            % RF spoiling
            rf.phaseOffset = mod(117*(iY^2+iY+2)*pi/180, 2*pi);
            adc.phaseOffset = rf.phaseOffset;
            
            % Excitation
            % Mark start of block group (= one TR) by adding label
            % (subsequent blocks in block group are not labelled).
            %seq.addBlock(rf, rfDelay);
            blockGroupID = 1;
            seq.addBlock(rf, mr.makeLabel('SET', 'LIN', blockGroupID));
            
            % Encoding
            seq.addBlock(mr.makeDelay(delayTE(c)));
            seq.addBlock(gxPre, ...
                mr.scaleGrad(gyPre, yStep), ...
                mr.scaleGrad(gzPre, zStep));
            if (iZ < 0)
                seq.addBlock(gx);
            else
                seq.addBlock(gx, adc);
            end

            % rephasing/spoiling
            seq.addBlock(gxSpoil, ...
                mr.scaleGrad(gyPre, -yStep), ...
                mr.scaleGrad(gzPre, -zStep));
            %seq.addBlock(gzSpoil);
            seq.addBlock(mr.makeDelay(delayTR(c)));
        end
    end
end
fprintf('\nSequence ready\n');

% Check sequence timing
[ok, error_report]=seq.checkTiming;
if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

% Visualise sequence and output for execution
Ndummy = length(TE)*Ny*nDummyZLoops;
seq.plot('TimeRange',[Ndummy+1 Ndummy+6]*TR(1), 'timedisp', 'ms')

seq.setDefinition('FOV', fov);
seq.setDefinition('Name', 'b0');
seq.write('b0.seq', false);

return

% visualize the 3D k-space (only makes sense for low-res, otherwise one sees nothing)
if Nx<=32
    tic;
    [kfa,ta,kf]=seq.calculateKspacePP();
    toc
    figure;plot3(kf(1,:),kf(2,:),kf(3,:));
    hold on;plot3(kfa(1,:),kfa(2,:),kfa(3,:),'r.');
end
