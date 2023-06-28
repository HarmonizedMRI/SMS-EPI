% SMS-EPI sequence in Pulseq
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

%% Paths
caipiPythonPath = '~/github/HarmonizedMRI/3DEPI/caipi/';

%% Define experimental parameters
sys = mr.opts('maxGrad', 50, 'gradUnit','mT/m', ...
              'maxSlew', 120, 'slewUnit', 'T/m/s', ...
              'rfDeadTime', 100e-6, ...
              'rfRingdownTime', 60e-6, ...
              'adcDeadTime', 40e-6, ...
              'adcRasterTime', 2e-6, ...
              'gradRasterTime', 10e-6, ...
              'blockDurationRaster', 10e-6, ...
              'B0', 3.0);

%timessi = 100e-6;    % start sequence interrupt (SSI) time (required delay at end of block group/TR)

voxelSize = [2.4 2.4 2.4]*1e-3;     % m
mb = 6;                             % multiband/SMS factor
Nx = 92; Ny = Nx; Nz = mb*10;       % Matrix size
fov = voxelSize .* [Nx Ny Nz];      % FOV (m)
TE = 30e-3;
alpha = 60;                         % flip angle (degrees)

nFrames = 10;                       % number of temporal frames (image volumes)
nDummyFrames = 0;                   % dummy frames to reach steady state

dwell = 4e-6;                       % ADC sample time (s). For GE, must be multiple of 2us.

% CAIPI sampling parameters
Ry = 1;                   % ky undersampling factor
Rz = mb;                  % kz undersampling factor
CaipiShiftZ = 2;
pf_ky = 0.7;             % partial Fourier factor along ky
etl = ceil(pf_ky*Ny/Ry);  % echo train length

nCyclesSpoil = 2;    % number of spoiler cycles, along x and z

rfTB  = 6;          % RF pulse time-bandwidth product
rfDur = 8e-3;       % RF pulse duration (s)

fatChemShift = 3.5*1e-6;                        % 3.5 ppm
fatOffresFreq = sys.gamma*sys.B0*fatChemShift;  % Hz

%% Excitation pulse
sysGE = toppe.systemspecs('maxGrad', sys.maxGrad/sys.gamma*100, ...   % G/cm
    'maxSlew', sys.maxSlew/sys.gamma/10, ...           % G/cm/ms
    'maxRF', 0.25);
sliceSep = fov(3)/mb;   % center-to-center separation between SMS slices (m)
[rf, gzRF, freq] = getsmspulse(alpha, voxelSize(3), rfTB, rfDur, ...
    mb, sliceSep, sysGE, sys, ...
    'doSim', false, ...    % Plot simulated SMS slice profile
    'type', 'st', ...     % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
    'ftype', 'ls');       % filter design. 'ls' = least squares

%% Get CAIPI sampling pattern (for one shot/echo train)
pyFile = [caipiPythonPath '/skippedcaipi_sampling.py'];
pyCmd = sprintf('python %s %d %d %d %d %d %d', ...
    pyFile, Ny, Nz, Ry, Rz, CaipiShiftZ, 1);
a = input('Press 1 to run Python script from Matlab, 2 to run offline:  ');
if a == 1
    system(pyCmd);
else
    fprintf('Open a terminal and run the following python command:\n\t%s\n', pyCmd);
    input('Then press Enter to continue');
end
load caipi

% kz and ky indeces (in units of deltak)
kzInds = double(indices((end-etl+1):end, 1));
kyInds = double(indices((end-etl+1):end, 2));

% kz encoding blip amplitude along echo train (multiples of deltak)
kzStep = diff(kzInds);
kyStep = diff(kyInds);


%% Define readout gradients and ADC event
deltak = 1./fov;

% start with the blips
gyBlip = mr.makeTrapezoid('y', sys, 'Area', max(abs(kyStep))*deltak(2)); 
gzBlip = mr.makeTrapezoid('z', sys, 'Area', max(abs(kzStep))*deltak(3)); 

% area and duration of the biggest blip
if gyBlip.area > gzBlip.area
    maxBlipArea = gyBlip.area;
    blipDuration = mr.calcDuration(gyBlip);
else
    maxBlipArea = gzBlip.area;
    blipDuration = mr.calcDuration(gzBlip);
end

% Readout trapezoid
% Limit gradient amplitude according to FOV and dwell time.
systmp = sys;
systmp.maxGrad = deltak(1)/dwell;
gro = mr.makeTrapezoid('x', systmp, 'Area', Nx*deltak(1) + maxBlipArea);

% ADC event
Tread = mr.calcDuration(gro) - blipDuration;
if mod(round(Tread*1e6)*1e-6, dwell)
    Tread = Tread - mod(Tread, dwell) + dwell;
end
adc = mr.makeAdc(round(Tread/dwell), sys, ...
    'Duration', Tread, ...
    'Delay', blipDuration/2);

% split blips at block boundary
[gyBlipUp, gyBlipDown] = mr.splitGradientAt(gyBlip, mr.calcDuration(gyBlip)/2);
gyBlipUp.delay = mr.calcDuration(gro) - mr.calcDuration(gyBlip)/2;
gyBlipDown.delay = 0;

[gzBlipUp, gzBlipDown] = mr.splitGradientAt(gzBlip, mr.calcDuration(gzBlip)/2);
gzBlipUp.delay = mr.calcDuration(gro) - mr.calcDuration(gzBlip)/2;
gzBlipDown.delay = 0;

% prephasers and spoilers
gxPre = mr.makeTrapezoid('x', sys, ...
    'Area', -gro.area/2);
Tpre = mr.calcDuration(gxPre);
gyPre = mr.makeTrapezoid('y', sys, ...
    'Area', (kyInds(1)-Ny/2+1/2)*deltak(2), ... 
    'Duration', Tpre);
gzPre = mr.makeTrapezoid('z', sys, ...
    'Area', -deltak(3), ...   % maximum PE2 gradient, max positive amplitude
    'Duration', Tpre);
gxSpoil = mr.makeTrapezoid('x', sys, ...
    'Area', Nx*deltak(1)*nCyclesSpoil);
gzSpoil = mr.makeTrapezoid('z', sys, ...
    'Area', Nx*deltak(1)*nCyclesSpoil);

%% Calculate delay to achieve desired TE
kyIndAtTE = find(kyInds-Ny/2 == min(abs(kyInds-Ny/2)));
minTE = mr.calcDuration(gzRF) - mr.calcDuration(rf)/2 - rf.delay + mr.calcDuration(gxPre) + ...
        (kyIndAtTE-0.5) * mr.calcDuration(gro);
TEdelay = floor((TE-minTE)/sys.blockDurationRaster) * sys.blockDurationRaster;

%% Assemble sequence
seq = mr.Sequence(sys);           

Nshots = Nz/mb;
kyStepMax = max(abs(kyStep));
kzStepMax = max(abs(kzStep));

        blockGroupID = 1;
        seq.addBlock(rf, gzRF, mr.makeLabel('SET', 'LIN', blockGroupID));
        if TE > minTE
            seq.addBlock(mr.makeDelay(TEdelay));
        end
        seq.addBlock(gxPre, gyPre, gzPre);
        seq.addBlock(gro, adc, ...
                     mr.scaleGrad(gyBlipUp, kyStep(1)/kyStepMax), ...
                     mr.scaleGrad(gzBlipUp, kzStep(1)/kzStepMax));
        for ie = 2:(etl-1)
            gybd = mr.scaleGrad(gyBlipDown, kyStep(ie-1)/kyStepMax);
            gybu = mr.scaleGrad(gyBlipUp, kyStep(ie)/kyStepMax);
            gybdu = mr.addGradients({gybd, gybu}, sys);
            gzbd = mr.scaleGrad(gzBlipDown, kzStep(ie-1)/kzStepMax);
            gzbu = mr.scaleGrad(gzBlipUp, kzStep(ie)/kzStepMax);
            gzbdu = mr.addGradients({gzbd, gzbu}, sys);
            seq.addBlock(adc, mr.scaleGrad(gro, (-1)^(ie-1)), gybdu, gzbdu);
        end
        seq.addBlock(adc, ...
                     mr.scaleGrad(gro, (-1)^(ie)), ...
                     mr.scaleGrad(gyBlipDown, kyStep(ie)/kyStepMax), ...
                     mr.scaleGrad(gzBlipDown, kzStep(ie)/kzStepMax));

% Check sequence timing
[ok, error_report]=seq.checkTiming;
if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% Inspect sequence

%seq.plot('blockrange', [1 20]);

% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D k-space plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
title('full k-space trajectory (k_x x k_y)');

return

% Calculate timing
%TEmin = rf.shape_dur/2 + rf.ringdownTime + mr.calcDuration(gxPre) ...
%      + adc.delay + Nx/2*dwell;
%delayTE = ceil((TE-TEmin)/seq.gradRasterTime)*seq.gradRasterTime;
%TRmin = mr.calcDuration(rf) + delayTE + mr.calcDuration(gxPre) ...
%      + mr.calcDuration(gx) + mr.calcDuration(gxSpoil);
%delayTR = ceil((TR-TRmin)/seq.gradRasterTime)*seq.gradRasterTime;

%% Loop over phase encodes and define sequence blocks
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

% Visualise sequence and output for execution
Ndummy = length(TE)*Ny*nDummyZLoops;
%seq.plot('TimeRange',[Ndummy+1 Ndummy+6]*TR(1), 'timedisp', 'ms')

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
