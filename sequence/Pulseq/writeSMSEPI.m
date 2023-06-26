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

%% Paths
caipiPythonPath = '~/github/HarmonizedMRI/3DEPI/caipi/';

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
TE = 30e-3;
alpha = 60;                         % flip angle (degrees)

dwell = 2e-6;                       % ADC sample time (s). For GE, must be multiple of 2us.

% CAIPI sampling parameters
Ry = 1;   % ky undersampling factor
Rz = mb;  % kz undersampling factor
CaipiShiftZ = 2;
pf_ky = 0.8;              % partial Fourier factor along ky
etl = ceil(pf_ky*Ny/Ry);  % echo train length

nCyclesSpoil = 2;    % number of spoiler cycles, along x and z

rfTB  = 6;          % RF pulse time-bandwidth product
rfDur = 8e-3;       % RF pulse duration (s)

fatChemShift = 3.5*1e-6;                        % 3.5 ppm
fatOffresFreq = sys.gamma*sys.B0*fatChemShift;  % Hz

Tpre = 0.5e-3;    % x/y/z prephasing gradient lob duration

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

%% Get CAIPI sampling pattern for each shot
pyFile = [caipiPythonPath '/skippedcaipi_sampling.py'];
pyCmd = sprintf('python %s %d %d %d %d %d %d', ...
    pyFile, Ny, Nz, Ry, Rz, CaipiShiftZ, 1);
a = 2; %input('Press 1 to run Python script from Matlab, 2 to run offline ');
if a == 1
    system(pyCmd);
else
    fprintf('Run the following python command:\n\t%s\n', pyCmd);
    input('Then press Enter to continue');
end
load caipi

% kz encoding blip amplitude along echo train (multiples of deltak)
Kzstep = diff(double(indices(1:etl,1) + 1));
Kystep = diff(double(indices(1:etl,2) + 1));


%% Define other gradients and ADC events

deltak = 1./fov;

% start with the blip
gyBlip = mr.makeTrapezoid('y', sys, 'Area', max(abs(Kystep))*deltak(2)); 
gzBlip = mr.makeTrapezoid('z', sys, 'Area', max(abs(Kzstep))*deltak(3)); 

% readout trapezoid and ADC (ramp sampling)
if gyBlip.area > gzBlip.area
    maxBlipArea = gyBlip.area;   % max blip size is along y
    blipDuration = mr.calcDuration(gyBlip);
else
    maxBlipArea = gzBlip.area;   % max blip size is along z
    blipDuration = mr.calcDuration(gzBlip);
end
gro = mr.makeTrapezoid('x', sys, 'Area', Nx*deltak(1) + maxBlipArea);
adc = mr.makeAdc(Nx, sys, ...
    'Duration', Nx*dwell, ...
    'Delay', blipDuration/2);

% prephasers and spoilers
gxPre = mr.makeTrapezoid('x', sys, ...
    'Area', -gro.area/2, ...
    'Duration', Tpre);
gyPre = mr.makeTrapezoid('y', sys, ...
    'Area', Ny*deltak(2)/2, ...   % maximum PE1 gradient, max positive amplitude
    'Duration', Tpre);
gzPre = mr.makeTrapezoid('z', sys, ...
    'Area', Nz*deltak(3)/2, ...   % maximum PE2 gradient, max positive amplitude
    'Duration', Tpre);
gxSpoil = mr.makeTrapezoid('x', sys, ...
    'Area', Nx*deltak(1)*nCyclesSpoil);
gzSpoil = mr.makeTrapezoid('z', sys, ...
    'Area', Nx*deltak(1)*nCyclesSpoil);

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
seq = mr.Sequence(sys);           
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
