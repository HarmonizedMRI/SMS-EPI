% Create 2D EPI fMRI scan in TOPPE (and later in Pulseq)
%
% Output:
%   TOPPE scan files:  modules.txt, scanloop.txt, tipdown.mod, and readout.mod
%
% TODO: create the corresponding Pulseq file

% Set paths to Pulseq and TOPPE libraries
%addpath ~/github/pulseq/matlab/         % +mr package
%addpath ~/github/toppeMRI/toppe/        % +toppe package
%addpath ~/github/toppeMRI/PulseGEq/     % +pulsegeq package (Pulseq <--> TOPPE conversion)
% paths for MZ
%addpath ~/pulseq_home/github/pulseq/matlab/
%addpath ~/pulseq_home/github/toppe/
%addpath ~/pulseq_home/github/PulseGEq/

%% Sequence parameters
ex.flip = 90;        % flip angle (degrees)
ex.thick = 0.4;      % slice thickness (cm)
ex.tbw = 8;          % time-bandwidth product
ex.dur = 4;          % msec
ex.nSpoilCycles = 8;   %  number of cycles of gradient spoiling across slice thickness
ex.spacing = 0.5;    % center-to-center slice separation (cm)

nslices = 20;
SLICES = [1:2:nslices 2:2:nslices];   % slice ordering (minimize slice crosstalk)

scandur = 30;        % seconds

fov = 22.4;          % cm
nx = 64; ny = 64;    % matrix size
espmin = 0.512;      % minimum echo spacing allow by scanner (ms) (on GE, see /usr/g/bin/epiesp.dat)

delay = 10;          % (ms) delay after RF pulse. Determines TE. 

gamma = 4.2576e3;   % Hz/Gauss

%% Set hardware limits
% NB! maxGrad must be equal to physical hardware limit since used to scale gradients.
% maxSlew can be less than physical limit (for reducing PNS)
system.ge = toppe.systemspecs('maxSlew', 20, 'slewUnit', 'Gauss/cm/ms', ...
	'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...  
	'maxRf', 0.25, 'rfUnit', 'Gauss');

%system.siemens = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
%    'MaxSlew', 150, 'SlewUnit', 'T/m/s', ...
%	 'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

%% Slice selective excitation (use TOPPE wrapper around John Pauly's SLR toolbox)
[ex.rf, ex.g, ex.freq] = toppe.utils.rf.makeslr(ex.flip, ex.thick, ex.tbw, ex.dur, ex.nSpoilCycles, ...
	'type', 'ex', ...   % 90 excitation. 'st' = small-tip
	'ftype', 'ls', ...  
	'ofname', 'tipdown.mod', ...
	'sliceOffset', ex.spacing, ...  % for calculating ex.freq
	'spoilDerate', 0.4, ...  % reduce slew by this much during spoiler/rewinder
	'system', system.ge);

%% EPI readout
res = fov/nx;          % spatial resolution (cm)
kmax = 1/(2*res);      % cycles/cm
area = kmax/gamma;     % G/cm * sec (area of each readout trapezoid)

% x/y prephaser
% reduce slew to reduce PNS
gpre = toppe.utils.trapwave2(area, system.ge.maxGrad, 0.8*system.ge.maxSlew/sqrt(2), system.ge.raster*1e3); % raster time in msec (sorry)
gpre = gpre(1:(end-1)); % remove 0 at end

% readout trapezoid
% Reduce maxGrad until desired echo spacing is reached.
% Allow ramp sampling, and violate Nyquist slightly near kmax.
for s = 1:-0.02:0.1
	mxg = s*system.ge.maxGrad;
	gx1 = toppe.utils.trapwave2(2*area, mxg, system.ge.maxSlew, system.ge.raster*1e3);
	if length(gx1)*system.ge.raster*1e3 > espmin
		break;
	end
end
gx1 = gx1(1:(end-1));  % remove 0 at end

% y blip
gyblip = toppe.utils.trapwave2(area/ny, system.ge.maxGrad, system.ge.maxSlew, system.ge.raster*1e3);

% gy waveform for 1st/last and other echoes
imax = find(gyblip == max(gyblip));
gyblipstart = gyblip(1:imax(1));  % first half of blip
gyblipend= gyblip(imax(2):end);
gy1 = [zeros(1,length(gx1)-length(gyblipstart)) gyblipstart]; % first echo
gyn = [gyblipend zeros(1,length(gx1)-length(gyblipend)-length(gyblipstart)) gyblipstart]; % other echoes
gylast = [gyblipend zeros(1,length(gx1)-length(gyblipend))]; % last echo

% put it all together and write to readout.mod
gx = [-gpre gx1];
gy = [-gpre gy1];
for iecho = 2:(ny-1)
	gx = [gx gx1*(-1)^(iecho+1)];
	gy = [gy gyn];
end
gx = [gx gx1*(-1)^(iecho+2) 0];  % add zero at end
gy = [gy gylast 0];
gx = toppe.utils.makeGElength(gx(:));
gy = toppe.utils.makeGElength(gy(:));
toppe.writemod('gx', gx, 'gy', gy, ...
	'system', system.ge, ...
	'ofname', 'readout.mod');

%% Create modules.txt 
% Entries are tab-separated
modFileText = ['' ...
'Total number of unique cores\n' ...
'2\n' ...
'fname	duration(us)	hasRF?	hasDAQ?\n' ...
'tipdown.mod	0	1	0\n' ...
'readout.mod	0	0	1' ];
fid = fopen('modules.txt', 'wt');
fprintf(fid, modFileText);
fclose(fid);

% determine sequence TR and number of frames
toppe.write2loop('setup', 'version', 3);
toppe.write2loop('tipdown.mod', 'textra', delay);
toppe.write2loop('readout.mod');
toppe.write2loop('finish');
trseq = toppe.getTRtime(1,2);       % sec
nframes = 2*ceil(scandur/(trseq*nslices)/2); % force to be even


%% Create scanloop.txt

% rf spoiling isn't really needed since T2 << TR but why not
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;

toppe.write2loop('setup', 'version', 3);   % Initialize scanloop.txt
for ifr = 1:nframes
	% turn off gy for odd/even echo calibration (first and 2nd-to-last frame)
	gyamp = 1.0 - any(ifr == [1 nframes-1]);

	% flip gx for 2D odd/even echo calibration (2nd and last frame)
	gxamp = (-1)^any(ifr == [2 nframes]);

	for isl = SLICES
		% excitation
	  	toppe.write2loop('tipdown.mod', 'RFamplitude', 1.0, ...
			'RFphase', rfphs, ...
			'textra', delay, ...
			'RFoffset', round((isl-0.5-nslices/2)*ex.freq) );  % Hz (slice selection)

	 	% readout
		% data is stored in 'slice', 'echo', and 'view' indeces. Will change to ScanArchive in future.
		toppe.write2loop('readout.mod', ...
			'Gamplitude', [gxamp gyamp 0]', ...
			'DAQphase', rfphs, ...
			'slice', isl, 'echo', 1, 'view', ifr, ...  
			'dabmode', 'on');

		% update rf phase (RF spoiling)
		rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
		rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
	end
end
toppe.write2loop('finish');

return;


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
