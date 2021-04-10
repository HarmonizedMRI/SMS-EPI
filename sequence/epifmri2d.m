% Create 2D EPI fMRI scan in TOPPE (and later in Pulseq)
%
% Later: create cross-vendor (GE and Siemens) SMS EPI fMRI scan 
%
% Outputs:
%  Pulseq scan file        smsepifmri.seq (for Siemens)
%  TOPPE scan files        smsepifmri.tgz (for GE)
%
% Usage:

% Set paths to Pulseq and TOPPE libraries
%addpath ~/github/pulseq/matlab/         % +mr package
%addpath ~/github/toppeMRI/toppe/        % +toppe package
%addpath ~/github/toppeMRI/PulseGEq/     % +pulsegeq package (Pulseq <--> TOPPE conversion)
% paths for MZ
%addpath ~/pulseq_home/github/pulseq/matlab/
%addpath ~/pulseq_home/github/toppe/
%addpath ~/pulseq_home/github/PulseGEq/

gamma = 4.2576e3;   % Hz/Gauss

%% Set hardware limits
% NB! maxGrad must be equal to physical hardware limit since used to scale gradients.
% maxSlew can be less than physical limit (for reducing PNS)
system.ge = toppe.systemspecs('maxSlew', 20, 'slewUnit', 'Gauss/cm/ms', ...
	'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...  
	'maxRf', 0.25, 'rfUnit', 'Gauss');

system.siemens = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
    'MaxSlew', 150, 'SlewUnit', 'T/m/s', ...
	 'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

%% Slice selective excitation (use TOPPE wrapper around John Pauly's SLR toolbox)
ex.flip = 90;        % degrees
ex.thick = 0.4;      % slice thickness (cm)
ex.spacing = 0.5;    % center-to-center slice separation (cm)
ex.tbw = 6;          % time-bandwidth product
ex.dur = 2;          % msec
nSpoilCycles = 8;   %  number of cycles of gradient spoiling across slice thickness
[ex.rf, ex.g, ex.freq] = toppe.utils.rf.makeslr(ex.flip, ex.thick, ex.tbw, ex.dur, nSpoilCycles, ...
	'type', 'ex', ...   % 90 excitation. 'st' = small-tip
	'ftype', 'ls', ...  
	'ofname', 'tipdown.mod', ...
	'sliceOffset', ex.spacing, ...  % for calculating ex.freq
	'system', system.ge);

%% EPI readout
fov=22.4;              % cm
nx=64; ny=64;          % matrix size
res = fov/nx;          % spatial resolution (cm)
kmax = 1/(2*res);      % cycles/cm
area = kmax/gamma;     % G/cm * sec (area of each readout trapezoid)
esp = 0.512;           % desired echo spacing (ms). Chosen so that 'forbidden' EPI spacings are not violated.

% x/y prephaser
gpre = toppe.utils.trapwave2(area, system.ge.maxGrad, system.ge.maxSlew, system.ge.raster*1e3); % raster time in msec (sorry)

% readout trapezoid
% reduce maxGrad until desired echo spacing is reached (allow ramp sampling)
for s = 1:-0.02:0.1
	mxg = s*system.ge.maxGrad;
	gx1 = toppe.utils.trapwave2(2*area, mxg, system.ge.maxSlew, system.ge.raster*1e3);
	if length(gx1)*system.ge.raster*1e3 > esp
		break;
	end
end

% y blip
gyblip = toppe.utils.trapwave2(area/ny, system.ge.maxGrad, system.ge.maxSlew, system.ge.raster*1e3);

% gy waveform for 1st and other echoes
imax = find(gyblip == max(gyblip));
gyblipstart = gyblip(1:imax(1));  % first half of blip
gyblipend= gyblip(imax(2):end);
gy1 = [zeros(1,length(gx1)-length(gyblipend)) gyblipstart]; % first echo
gyn = [gyblipend zeros(1,length(gx1)-length(gyblipend)-length(gyblipstart)) gyblipstart]; % other echoes

% put it all together
gx = [-gpre gx1];
gy = [-gpre gy1];
for iecho = 2:ny
	gx = [gx gx1*(-1)^(iecho+1)];
	gy = [gy gy1];
end

% reshape
g = [gx(:) gy(:)];


return

% Phase blip in shortest possible time
blip_dur = ceil(2*sqrt(deltak/system.siemens.maxSlew)/10e-6/2)*10e-6*2; % we round-up the duration to 2x the gradient raster time
% the split code below fails if this really makes a trpezoid instead of a triangle...
gy = mr.makeTrapezoid('y',system.siemens,'Area',deltak,'Duration',blip_dur);

return;

% Create modules.txt (do it here since this file is used for intermediate steps below)
% Entries are tab-separated
modFileText = ['' ...
'Total number of unique cores\n' ...
'2\n' ...
'fname	duration(us)	hasRF?	hasDAQ?\n' ...
'readout.mod	0	0	1\n' ...
'tipdown.mod	0	1	0' ];
fid = fopen('modules.txt', 'wt');
fprintf(fid, modFileText);
fclose(fid);

% Create a new Pulseq sequence object
%seq = mr.Sequence(siemens.system);         

% Design rf waveform
nSpoilCycles = 1e-3;   % just has to be small enough so that no spoiler is added at beginning of waveform in makeslr()
[ex.rf, ex.g] = toppe.utils.rf.makeslr(ex.flip, ex.thick, ex.tbw, ex.dur, nSpoilCycles, ...
	'type', 'ex', ...   % 90 excitation. 'st' = small-tip
	'ftype', 'ls', ...  
	'system', arg.system, ...
	'writeModFile', false);

% Design readout gradients
% Remember: makegre() creates y phase-encodes based on isotropic in-plane resolution,
% so need to scale gy in scanloop.txt accordingly (see below).
dzDummy = 1e3;   % z voxel size. Just needs to be large enough so the trapezoid doesn't extend overall waveform duration.
[acq.gx, acq.gy, acq.gz] = toppe.utils.makegre(fov, nx, dzDummy, ...
	'system', limits.design, 'ncycles', nCycleSpoil, ...
	'ofname', 'tmp.mod', ...
	'slewDerate', 0.7, ...   % reduce PNS
	'oprbw', oprbw);

% Create waveforms suitable for Pulseq by converting units and interpolating
% discard zeros at beginning and end of rf waveform
I = find(abs(ex.rf) == 0);
iStart = find(diff(I)>1) + 1;
iStop = I(iStart+1);
siemens.ex.rf = pulsegeq.rf2pulseq(ex.rf(iStart:iStop), ge.system.raster, seq);  % Gauss -> Hz; 4us -> 1us.
siemens.ex.rfdelay = roundToRaster(iStart * ge.system.raster, siemens.system.gradRasterTime);
siemens.ex.gdelay = max(0, siemens.system.rfDeadTime - siemens.ex.rfdelay);
siemens.ex.g  = pulsegeq.g2pulseq(ex.g, ge.system.raster, seq);    % Gauss/cm -> Hz/m; 4us -> 10us
siemens.acq.gx = pulsegeq.g2pulseq(acq.gx, ge.system.raster, seq); % readout gradient
siemens.acq.gy = pulsegeq.g2pulseq(acq.gy, ge.system.raster, seq); % phase-encode gradient
siemens.acq.gz = pulsegeq.g2pulseq(acq.gz, ge.system.raster, seq); % partition-encode gradient

% Write .mod files for TOPPE
% For readout, load 'tmp.mod' created above, swap gx/gz, and write to 'readout.mod'
toppe.writemod('rf', toppe.utils.makeGElength(ex.rf), ...
	'gz', toppe.utils.makeGElength(ex.g), ...
	'ofname', 'tipdown.mod', ...
	'system', ge.system);

[rf,gx,gy,gz,desc,paramsint16,paramsfloat] = toppe.readmod('tmp.mod');
system('rm tmp.mod');
toppe.writemod('gx', gz, 'gy', gy, 'gz', gx, 'ofname', 'readout.mod', ...
	'system', ge.system, ...
	'desc', desc, 'hdrints', paramsint16, 'hdrfloats', oprbw);

% Determine delays to achieve desired TR.
% For TOPPE we do this by writing a short scanloop.txt file and calculating the resulting minimum TR from that.
toppe.write2loop('setup', 'version', 3);
toppe.write2loop('tipdown.mod');
toppe.write2loop('readout.mod');
toppe.write2loop('finish');
tr_min = toppe.getTRtime(1,2)*1e3;       % msec
ge.delay = TR - tr_min;                  % msec

tmp.rf = mr.makeArbitraryRf(siemens.ex.rf, ex.flip/180*pi, 'system', siemens.system, 'delay', siemens.ex.rfdelay);
tmp.readout = mr.makeArbitraryGrad('z', siemens.acq.gx, siemens.system);
tr_min = mr.calcDuration(tmp.rf) + mr.calcDuration(tmp.readout);   % sec
siemens.delay = roundToRaster(TR*1e-3-tr_min, siemens.system.gradRasterTime);       % sec

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
pulseq.acq.gz = mr.makeArbitraryGrad('z', siemens.acq.gx, siemens.system); 
pulseq.ex.rf = mr.makeArbitraryRf(siemens.ex.rf, ex.flip/180*pi, ...
	'system', siemens.system, 'delay', siemens.ex.rfdelay);
pulseq.ex.g  = mr.makeArbitraryGrad('z', siemens.ex.g, siemens.system, 'delay', siemens.ex.gdelay);

% Scan loop
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;

toppe.write2loop('setup', 'version', 3);   % Initialize scanloop.txt
for iz = 1:nz         
	for iy = (-nDisdaq+1):ny 
		if iy > 0   
			dabmode = 'on';
		else
			dabmode = 'off';
		end

		% y/z phase-encode amplitudes, scaled to (-1,1)
		a_gy = (iy > 0) * ((iy-1+0.5)-ny/2)/(ny/2) * ny/nx;  % note factor of ny/nx
		a_gz = (iy > 0) * ((iz-1+0.5)-nz/2)/(nz/2);

		% rf excitation (TOPPE)
	  	toppe.write2loop('tipdown.mod', 'RFamplitude', 1.0, 'RFphase', rfphs);

		% rf excitation (Pulseq)
		pulseq.ex.rf.phaseOffset = rfphs;
		seq.addBlock(pulseq.ex.rf, pulseq.ex.g);

	 	% readout (TOPPE)
		toppe.write2loop('readout.mod', ...
			'Gamplitude', [a_gz a_gy 1]', ...
			'DAQphase', rfphs, ...
			'slice', iz, 'echo', 1, 'view', max(iy,1), ...  % data is stored in 'slice', 'echo', and 'view' indeces. Will change to ScanArchive in future.
			'textra', ge.delay, ...                         % dead time at end of module (msec)
			'dabmode', dabmode);

		% readout (Pulseq)
		gx = mr.makeArbitraryGrad('x', a_gz*siemens.acq.gz, siemens.system); 
		gy = mr.makeArbitraryGrad('y', a_gy*siemens.acq.gy, siemens.system); 
		pulseq.adc.phaseOffset = rfphs;   % radians
		seq.addBlock(gx, gy, pulseq.acq.gz, pulseq.adc); 
		seq.addBlock(mr.makeDelay(siemens.delay));  % adding delay event to readout block does not produce desired TR. Is this the intended Pulseq behavior?

		% update rf phase (RF spoiling)
		rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
		rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
	end
end
toppe.write2loop('finish');

% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

% Write to Pulseq file
seq.setDefinition('FOV', [fov fov ex.thick]*1e-2);   % m
seq.setDefinition('Name', '2D MIP of SMS slice profile');
seq.write('SMSprofile.seq');
% parsemr('SMSprofile.seq');

% Write TOPPE files to a tar file
system('tar czf SMSprofile.tgz modules.txt scanloop.txt tipdown.mod readout.mod');

% Display (part of) sequences
nModsPerTR = 2;    % number of TOPPE modules per TR
nTR = ny/2;        % number of TRs to display
nStart = nModsPerTR * floor(nDisdaq+ny/2-nTR/2);
toppe.plotseq(nStart, nStart + nTR*nModsPerTR);
tStart = nStart/nModsPerTR*TR*1e-3;    % sec
tStop = tStart + nTR*TR*1e-3;          % sec
seq.plot('timeRange', [tStart tStop]);

% Display TOPPE sequence in loop/movie mode
fprintf('Displaying sequence...');
%figure; toppe.playseq(nModsPerTR, 'drawpause', false);  
fprintf('done\n');
