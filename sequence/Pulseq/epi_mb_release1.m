% using the demo echo train, change the blip_dur so far
% rf, refscan all in one code
% caipi using the one dimension
% ref using the full scan of the EPI

%% settings
sys = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
    'MaxSlew', 150, 'SlewUnit', 'T/m/s', ... 
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

% basic parameters
seq=mr.Sequence(sys);           % Create a new sequence object
fov=256e-3; Nx=128; Ny=Nx;      % Define FOV and resolution
alpha=90;                       % flip angle
slicethickness=3e-3;            % slice
%TR=21e-3;                      % ignore TR, go as fast as possible
%TE=60e-3;                      % ignore TE, go as fast as possible
% TODO: run this EPI sequence with and without PE, to calibrate the delay
% TODO: change MaxGrad/MaxSlew/roDuration to see what happends to the
% stimulation (e.g. 80/200/500)

% more in-depth parameters
pe_enable=1;                    % a flag to quickly disable phase encoding (1/0) as needed for the delay calibration
rfDuration=3e-3;
roDuration=640e-6;              % not all values are possible, watch out for the checkTiming output
%% blocks: rf, blips, readout, spoiler
% Create alpha-degree slice selection pulse and corresponding gradients 
[rf, gz, gzReph] = mr.makeSincPulse(alpha*pi/180,'Duration',rfDuration,...
    'SliceThickness',slicethickness,'apodization',0.42,'timeBwProduct',4,'system',sys);

% Define other gradients and ADC events
deltak=1/fov; % Pulseq default units for k-space are inverse meters
kWidth = Nx*deltak;

% start with the blip
blip_dur = ceil(2*sqrt(deltak/sys.maxSlew)/sys.gradRasterTime/2)*sys.gradRasterTime*2; % we round-up the duration to 2x the gradient raster time
blip_dur = blip_dur*2;
blip_dur = 160e-6;
gyBlip = mr.makeTrapezoid('y',sys,'Area',-deltak,'Duration',blip_dur); % we use negative blips to save one k-space line on our way towards the k-space center

% readout gradient is a truncated trapezoid with dead times at the beginnig
% and at the end, each equal to a half of blip_dur
% the area between the blips should be equal to kWidth
% we do a two-step calculation: we first increase the area assuming maximum
% slew rate and then scale down the amlitude to fix the area 
extra_area=blip_dur/2*blip_dur/2*sys.maxSlew;
gx = mr.makeTrapezoid('x',sys,'Area',kWidth+extra_area,'duration',roDuration+blip_dur);
actual_area=gx.area-gx.amplitude/gx.riseTime*blip_dur/2*blip_dur/2/2-gx.amplitude/gx.fallTime*blip_dur/2*blip_dur/2/2;
gx=mr.scaleGrad(gx,kWidth/actual_area);
%adc = mr.makeAdc(Nx,'Duration',roDurtion,'Delay',blip_dur/2,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'system',sys); % if no 'Duration' is provided shortest possible duration will be used
gyPre = mr.makeTrapezoid('y','Area',(Ny/2-1)*deltak,'system',sys);

% calculate ADC - it is quite trickly
% we use ramp sampling, so we have to calculate the dwell time and the
% number of samples, which are will be quite different from Nx and
% readoutTime/Nx, respectively. 
adcDwellNyquist=deltak/gx.amplitude; % dwell time on the top of the plato
% round-down dwell time to sys.adcRasterTime (100 ns)
adcDwell=floor(0.5*adcDwellNyquist/sys.adcRasterTime)*sys.adcRasterTime;
adcSamples=floor(roDuration/adcDwell/4)*4; % on Siemens the number of ADC samples need to be divisible by 4
adc = mr.makeAdc(adcSamples,'Dwell',adcDwell);
% realign the ADC with respect to the gradient
time_to_center=adc.dwell*((adcSamples-1)/2+0.5); % Pulseq (andsiemens) define the samples to happen in the center of the dwell period
adc.delay=round((gx.riseTime+gx.flatTime/2-time_to_center)/sys.rfRasterTime)*sys.rfRasterTime; 
          % above we adjust the delay to align the trajectory with the gradient.
          % We have to aligh the delay to seq.rfRasterTime (1us) 
          % this rounding actually makes the sampling points on odd and even readouts
          % to appear misalligned. However, on the real hardware this misalignment is
          % much stronger anyways due to the grdient delays

% finish the blip gradient calculation
% split the blip into two halves and produce a combined synthetic gradient
gyBlip_parts = mr.splitGradientAt(gyBlip, blip_dur/2, sys);
[gyBlip_up,gyBlip_down,~]=mr.align('right',gyBlip_parts(1),'left',gyBlip_parts(2),gx);
% now for inner echos create a special gy gradient, that will ramp down to 0, stay at 0 for a while and ramp up again
gyBlip_down_up=mr.addGradients({gyBlip_down, gyBlip_up}, sys);

% pe_enable support
gyBlip_up=mr.scaleGrad(gyBlip_up,pe_enable);
gyBlip_down=mr.scaleGrad(gyBlip_down,pe_enable);
gyBlip_down_up=mr.scaleGrad(gyBlip_down_up,pe_enable);
gyPre=mr.scaleGrad(gyPre,pe_enable);

% gradient spoiling
gzSpoil=mr.makeTrapezoid('z','Area',4/slicethickness,'system',sys); % 4 cycles over the slice thickness

% skip timing (TE/TR calculation), we'll accept the shortest TE/TR
%% blocks to caipi

%% multiband settings
Bandsep = [0 -7 -14 7 14 21];
nBand = 0;

% trials and errors, blip duration is too small
%rf_6 = mrz.makeMBPulse_less(rf.signal,'n_bands',0,'timeBwProduct',4,'band_sep',[0 -7 -14 7 14 21]*3e-3/slicethickness);%phs is for no centered slice
%rf_61 = mrz.makeMBPulse(rf.signal,'n_bands',6,'timeBwProduct',4,'band_sep',7*3e-3/slicethickness,'phs_0_pt','phs_mod');


rf1 = rf;% rf1 for multiband
%rf1.signal = rf_6;

if nBand == 0
    nBand = length(Bandsep);
    kz = 1/abs(diff(Bandsep(1:2))*slicethickness*nBand);%evenly distribute nBand
else
    kz = 1/(abs(Bandsep)*slicethickness*nBand);
end

% max Area
% 1/(170*sys.gamma*0.25*blip_dur^2)%compared with 1/(nBands*thick*band_sep)

gzBlip = mr.makeTrapezoid('z',sys,'Area',kz,'Duration',blip_dur,'maxGrad',38*sys.gamma,'maxSlew',170*sys.gamma); 

gzBlip2 = mr.makeTrapezoid('z',sys,'Area',kz*(nBand-1),'Duration',blip_dur,'maxGrad',38*sys.gamma,'maxSlew',170*sys.gamma); 
        %% undersampling strategy 1
        %could be negeralized to the strategy 3
gzBlip_parts = mr.splitGradientAt(gzBlip, blip_dur/2, sys);
[gzBlip_l,gzBlip_r,~]=mr.align('right',gzBlip_parts(1),'left',gzBlip_parts(2),gx);
% \_______
%  ¯¯¯¯¯¯¯\
gzBlip_neg=mr.addGradients({gzBlip_r, mr.scaleGrad(gzBlip_l,-1)}, sys);%figure,plot(gzBlip_neg.tt,gzBlip_neg.waveform)
%  _______/
% /¯¯¯¯¯¯¯
gzBlip_pos=mr.addGradients({mr.scaleGrad(gzBlip_r,-1), gzBlip_l}, sys);
        %% undersampling strategy caipi 3
        % [0 2*pi/3 -2*pi/3]
        % [0 pi/3 2/3*pi 3/3pi 4/3pi 5/3pi]
gzBlip_parts1 = mr.splitGradientAt(gzBlip, blip_dur/2, sys);
gzBlip_parts2 = mr.splitGradientAt(gzBlip2, blip_dur/2, sys);
[gzBlip_l1,gzBlip_r1,~]=mr.align('right',gzBlip_parts1(1),'left',gzBlip_parts1(2),gx);
[gzBlip_l2,gzBlip_r2,~]=mr.align('right',gzBlip_parts2(1),'left',gzBlip_parts2(2),gx);%figure,plot(gzBlip_l2.tt,gzBlip_l2.waveform)
% /¯¯¯¯¯¯¯\
gzBlip_neg1=mr.addGradients({mr.scaleGrad(gzBlip_r1,-1), mr.scaleGrad(gzBlip_l1,-1)}, sys);%figure,plot(gzBlip_neg1.tt,gzBlip_neg1.waveform)
%  _______/
% /¯¯¯¯¯¯¯
gzBlip_pos1=mr.addGradients({mr.scaleGrad(gzBlip_r1,-1), mr.scaleGrad(gzBlip_l2,1)}, sys); 
% \_______
%  ¯¯¯¯¯¯¯\
gzBlip_neg2=mr.addGradients({mr.scaleGrad(gzBlip_r2,1), mr.scaleGrad(gzBlip_l1,-1)}, sys);


%% assembly blocks
% #1, reference scan (full scan, selective slice index scan)
% #2, autocalibrate (in the sequence)
% #3, sms

% so far, too many blocks, not combined as one yet
seq.addBlock(mr.makeDelay(0.2));
refs = [-14 -7 0 7 14 21];%slice separation
    %% #1 reference scan
% rf.freqOffset = -gz.amplitude*3e-3*18;
% rf.phaseOffset = 2;
% define sequence blocks
for iS = 1:length(refs)
    rf.freqOffset = -gz.amplitude*slicethickness*refs(iS);
    seq.addBlock(rf,gz);
    seq.addBlock(mr.align('left',gyPre,gzReph,'right',gxPre));
    for i=1:Ny % loop over phase encodes
        if i==1
            seq.addBlock(gx,gyBlip_up,adc); % Read the first line of k-space with a single half-blip at the end
        elseif i==Ny
            seq.addBlock(gx,gyBlip_down,adc); % Read the last line of k-space with a single half-blip at the beginning
        else
            seq.addBlock(gx,gyBlip_down_up,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
        end 
        gx = mr.scaleGrad(gx,-1);   % Reverse polarity of read gradient),'Duration',mr.calcDuration(gxPre),'system',sys);
    end
    seq.addBlock(gzSpoil);
     seq.addBlock(mr.makeDelay(0.05));%for Trio, larger than 4 slices,
%      should add this...  for Nx Ny 128
end
    %% #3 normal SMS

seq.addBlock(rf1,gz);
seq.addBlock(mr.align('left',gyPre,gzReph,'right',gxPre));
for i=1:Ny % loop over phase encodes
    if i==1
        seq.addBlock(gx,gyBlip_up,adc); % Read the first line of k-space with a single half-blip at the end
    elseif i==Ny
        seq.addBlock(gx,gyBlip_down,adc); % Read the last line of k-space with a single half-blip at the beginning
    else
        seq.addBlock(gx,gyBlip_down_up,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
    end 
    gx = mr.scaleGrad(gx,-1);   % Reverse polarity of read gradient),'Duration',mr.calcDuration(gxPre),'system',sys);
end
seq.addBlock(gzSpoil);
   
    %% #4 SMS caipi strategy 1
    %  X  X
    %   X  X
    %    X  X
% Caipirinha
% seq.addBlock(mr.makeDelay(0.1));
seq.addBlock(rf1,gz);
seq.addBlock(mr.align('left',gyPre,gzReph,'right',gxPre));
for i=1:Ny % loop over phase encodes
    if i==1
        step = 0;
        seq.addBlock(gx,gyBlip_up,mr.scaleGrad(gzBlip_l1,-1),adc); % Read the first line of k-space with a single half-blip at the end
    elseif i==Ny
        step = step + 1;
        %check the gz moment
        switch mod(step,nBand)
            case num2cell(1:nBand-2)
                seq.addBlock(gx,gyBlip_down,mr.scaleGrad(gzBlip_r1,-1),adc);
            case nBand-1
                seq.addBlock(gx,gyBlip_down,mr.scaleGrad(gzBlip_r1,-1),adc);
            case 0
                seq.addBlock(gx,gyBlip_down,mr.scaleGrad(gzBlip_r2,1),adc);
        end
    else
        step = step+1;
    switch mod(step,nBand)
    case num2cell(1:nBand-2)
        seq.addBlock(gx,gyBlip_down_up,gzBlip_neg1,adc);
    case nBand-1
        seq.addBlock(gx,gyBlip_down_up,gzBlip_pos1,adc);
    case 0
        seq.addBlock(gx,gyBlip_down_up,gzBlip_neg2,adc);
    end

    end 
    gx = mr.scaleGrad(gx,-1);   % Reverse polarity of read gradient),'Duration',mr.calcDuration(gxPre),'system',sys);
end
seq.addBlock(gzSpoil);


%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% prepare sequence export
seq.setDefinition('FOV', [fov fov slicethickness]);
seq.setDefinition('Name', 'epi_2mb');
%% write
seq.write('epi_7slc6_caipi.seq')       % Write to pulseq file
%seq.install('siemens');

%% plot sequence and k-space diagrams

%seq.plot();
seq.plot('timeDisp','us','showBlocks',1); %detailed view

% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D k-space plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
title('full k-space trajectory (k_x x k_y)');

%% PNS calc

[pns_ok, pns_n, pns_c, tpns]=seq.calcPNS('idea/asc/MP_GPA_K2309_2250V_951A_AS82.asc'); % prisma
%[pns_ok, pns_n, pns_c, tpns]=seq.calcPNS('idea/asc/MP_GPA_K2309_2250V_951A_GC98SQ.asc'); % aera-xq
%[pns_ok, pns_n, pns_c, tpns]=seq.calcPNS('idea/asc/MP_GPA_K2298_2250V_793A_SC72CD_EGA.asc'); % TERRA-XR

if (pns_ok)
    fprintf('PNS check passed successfully\n');
else
    fprintf('PNS check failed! The sequence will probably be stopped by the Gradient Watchdog\n');
end

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  
rep = seq.testReport;
fprintf([rep{:}]);
