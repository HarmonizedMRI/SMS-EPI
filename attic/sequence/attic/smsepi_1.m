% Create SMS EPI fMRI scan in TOPPE (and later in Pulseq)
%
% Output:
%   TOPPE scan files:  modules.txt, scanloop.txt, tipdown.mod, and readout.mod
%
% TODO: create the corresponding Pulseq file

%% Sequence parameters

% parameters common to all sequences in this folder (SMS/2D EPI, 3D GRE coil calibration scan)
[seq, sys] = getparams;
fov = seq.fov;              % in-plane field of view (cm)
nx = seq.nx; ny = seq.ny;   % matrix size
gamma = sys.ge.gamma;    % Hz/Gauss
dt = sys.ge.raster;      % sec

% SMS excitation
ex.flip = 30;        % flip angle (degrees)
ex.type = 'st';      % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
ex.ftype = 'ls';     
ex.tbw = 6;          % time-bandwidth product
ex.dur = 4;          % msec
ex.sliceSep = seq.slThick*6;     % center-to-center separation between SMS slices (cm)
mb = 6;          % sms/multiband factor (number of simultaneous slices)

% timing
delay.postrf = 5.0;    % (ms) delay after RF pulse. Determines TE. 
scandur = 1*20;       % seconds
tr = 300;             % (ms) If tr < minimum seq tr, minimum tr is calculated and used

% slice ordering
nshots = round(ex.sliceSep/seq.slThick);    % number of TRs for contiguous slice coverage
SHOTS = [1:2:nshots 2:2:nshots];    % slice ordering (minimize slice crosstalk)

nslices = mb*nshots;
zfov = seq.slThick*nslices;

fbesp = sys.ge.forbiddenEspRange;   % ms

%% SMS excitation waveforms
tmp = sys.ge;
tmp.maxSlew = 8;   % G/cm/ms. Reduce PNS during slice rephaser.
[ex.rf, ex.g, freq] = getsmspulse(ex.flip, seq.slThick, ex.tbw, ex.dur, mb, ex.sliceSep, tmp, ...
	'doSim', true, ...   % Plot simulated SMS slice profile
	'type', ex.type, ...
	'ftype', ex.ftype);

% write TOPPE module
toppe.writemod(sys.ge, 'rf', ex.rf, 'gz', ex.g, ...
	'ofname', 'tipdown.mod' );

freq = freq/ex.sliceSep*seq.slThick; % frequency offset for z shift of seq.slThick 

%% EPI readout waveforms
mxs = 0.55*sys.ge.maxSlew;
[gx, gy, gz, esp] = getepireadout(fov, nx, ny, ...
	sys.ge.maxGrad, mxs, dt*1e3, fbesp, ...
	mb, ex.sliceSep);

% write TOPPE module
toppe.writemod(sys.ge, 'gx', gx, 'gy', gy, 'gz', gz, ...
	'ofname', 'readout.mod' );

%% Gradient spoiler waveform
mxs = 7;  % Gauss/cm/ms. Lower to reduce PNS.
mxg = sys.ge.maxGrad;  % Gauss/cm
gcrush = toppe.utils.makecrusher(seq.nSpoilCycles, seq.slThick, sys.ge, 0, mxs, mxg);

% write TOPPE module
toppe.writemod(sys.ge, 'gx', gcrush, 'gy', gcrush, 'ofname', 'spoil.mod');

%% Create modules.txt for TOPPE 
% Entries are tab-separated
modFileText = ['' ...
'Total number of unique cores\n' ...
'3\n' ...
'fname	duration(us)	hasRF?	hasDAQ?\n' ...
'tipdown.mod	0	1	0\n' ...
'readout.mod	0	0	1\n' ...
'spoil.mod	0	0	0' ];
fid = fopen('modules.txt', 'wt');
fprintf(fid, modFileText);
fclose(fid);

%% Calculate delay to achieve desired TR
toppe.write2loop('setup', sys.ge, 'version', 3);
toppe.write2loop('tipdown.mod', sys.ge, 'textra', delay.postrf);
toppe.write2loop('readout.mod', sys.ge);
toppe.write2loop('spoil.mod', sys.ge);
toppe.write2loop('finish', sys.ge);
trseq = toppe.getTRtime(1,2,sys.ge)*1e3;    % sequence TR (ms)
if tr < trseq
	tr = trseq;
	fprintf('Using minimum tr (%.1f ms)\n', round(trseq));
end
	
delay.postreadout = tr-trseq;  % (ms) delay after readout

% number of frames
trvol = tr*nshots;  % ms
nframes = 2*ceil(scandur*1e3/trvol/2);  % force to be even

%% Prepare Pulseq sequence object and waveforms

if 0
% Create a new Pulseq sequence object
seq = mr.Sequence(siemens.system);         

% Create waveforms suitable for Pulseq by converting units and interpolating.
% Discard zeros at beginning and end of rf waveform.
I = find(abs(ex.rf) == 0);
iStart = find(diff(I)>1) + 1;
iStop = I(iStart+1);
siemens.ex.rf = pulsegeq.rf2pulseq(ex.rf(iStart:iStop), ge.system.raster, seq);  % Gauss -> Hz; 4us -> 1us.
siemens.ex.rfdelay = roundtoraster(iStart * ge.system.raster, siemens.system.gradRasterTime);
siemens.ex.gdelay = max(0, siemens.system.rfDeadTime - siemens.ex.rfdelay);
siemens.ex.g  = pulsegeq.g2pulseq(ex.g, ge.system.raster, seq);    % Gauss/cm -> Hz/m; 4us -> 10us
siemens.acq.gx = pulsegeq.g2pulseq(acq.gx, ge.system.raster, seq); % readout gradient
siemens.acq.gy = pulsegeq.g2pulseq(acq.gy, ge.system.raster, seq); % phase-encode gradient
siemens.acq.gz = pulsegeq.g2pulseq(acq.gz, ge.system.raster, seq); % partition-encode gradient

tmp.rf = mr.makeArbitraryRf(siemens.ex.rf, ex.flip/180*pi, 'system', siemens.system, 'delay', siemens.ex.rfdelay);
tmp.readout = mr.makeArbitraryGrad('z', siemens.acq.gx, siemens.system);
tr_min = mr.calcDuration(tmp.rf) + mr.calcDuration(tmp.readout);   % sec
siemens.delay = roundtoraster(TR*1e-3-tr_min, siemens.system.gradRasterTime);       % sec

% Create Pulseq adc object
[~,~,~,~,~,paramsint16] = toppe.readmod('readout.mod');
nPre = paramsint16(1);  % number of samples in gradient pre-winder and ramp to plateau
nPlateau = paramsint16(2);
acq.preDelay = nPre*ge.system.raster;            % sec
acq.flat = nPlateau*ge.system.raster;            % duration of flat portion of readout (sec)
siemens.acq.N = 2*round(acq.flat/siemens.system.gradRasterTime/2);   % number of readout samples (Siemens)
siemens.acq.dur = siemens.acq.N*siemens.system.gradRasterTime;
pulseq.adc = mr.makeAdc(siemens.acq.N, 'Duration', siemens.acq.dur, 'Delay', acq.preDelay, 'system', siemens.system);

% Create other Pulseq objects that don't need updating in scan loop (except phase)
pulseq.acq.gx = mr.makeArbitraryGrad('x', siemens.acq.gx, siemens.system); 
pulseq.ex.rf = mr.makeArbitraryRf(siemens.ex.rf, ex.flip/180*pi, ...
	'system', siemens.system, 'delay', siemens.ex.rfdelay);
pulseq.ex.g  = mr.makeArbitraryGrad('z', siemens.ex.g, siemens.system, 'delay', siemens.ex.gdelay);
end


%% Create scanloop.txt
% rf spoiling isn't really needed since T2 << TR but why not
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;

toppe.write2loop('setup', sys.ge, 'version', 4);   % Initialize scanloop.txt
for ifr = 1:nframes
	for ish = SHOTS
		% excitation
		f = round((ish-0.5-nshots/2)*freq);  % frequency offset (Hz) for slice shift
	  	toppe.write2loop('tipdown.mod', sys.ge, ...
            'RFamplitude', 1.0, ...
			'RFphase', rfphs, ...
			'textra', delay.postrf, ...
			'RFoffset', f );  % shift all slices

	 	% readout
		% data is stored in 'slice', 'echo', and 'view' indeces. Will change to ScanArchive in future.
		toppe.write2loop('readout.mod', sys.ge, ...
			'Gamplitude', [1 1 1]', ...
			'DAQphase', rfphs, ...
			'slice', ish, 'echo', 1, 'view', ifr, ...  
			'dabmode', 'on');
		
		% spoiler
		% set sign of gx/gy so they add to readout gradient first moment
		toppe.write2loop('spoil.mod', sys.ge, ...
			'textra', delay.postreadout, ...
			'Gamplitude', [sign(sum(gx)) sign(sum(gy)) 0]');

		% update rf phase (RF spoiling)
		rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
		rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
	end
end
toppe.write2loop('finish', sys.ge);

numModulesPerTR = 3;
figure; toppe.plotseq(1, numModulesPerTR, sys.ge, 'drawpause', false);

tar('smsepi.tar', {'*.txt', '*.mod', '*.m'});

fprintf('Created SMS EPI TOPPE files (smsepi.tar)\n');
fprintf('Matrix: [%d %d %d]; %.2f cm iso resolution; FOV: [%.1f %.1f %.1f] cm\n', ...
	nx, ny, nslices, fov/nx, fov, fov, zfov);
fprintf('SMS factor: %d; SMS slice separation: %.2f cm\n', mb, ex.sliceSep);
fprintf('Sequence TR: %.1f ms; Volume (frame) TR: %.1f ms\n', tr, trvol);
fprintf('Echo spacing: %.3f ms\n', esp);
fprintf('Number of frames: %d\n', nframes);
fprintf('To set TE, change delay.postrf in this script (currently set to %.1f ms)\n', delay.postrf);
fprintf('To plot first TR:       >> toppe.plotseq(1, 3, sys.ge, ''drawpause'', false);\n');
fprintf('To display scan loop:     >> toppe.playseq(3, sys.ge);       \n');
fprintf('To plot all .mod files: >> toppe.plotmod(''all''); \n');

return;


% simulate and plot slice profile
m0 = [0 0 1];        % initial magnetization
z = linspace(-fov/2, fov/2, 500);   % spatial locations (cm)
T1 = 1000;           % ms
T2 = 100;            % ms
dt = 4e-3;           % ms
figure;
[rf,gx,gy,gz] = toppe.readmod('tipdown.mod');
[m] = toppe.utils.rf.slicesim([0 0 1], rf, gz, dt, z, T1, T2, true);
subplot(132); title('simulated slice profile');

return;

