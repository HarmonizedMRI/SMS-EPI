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
ex.flip = 45;        % flip angle (degrees)
ex.type = 'st';      % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
ex.thick = 0.3;      % slice thickness (cm)
ex.spacing = 1.0;    % center-to-center slice separation (cm)
ex.tbw = 8;          % time-bandwidth product
ex.dur = 4;          % msec
ex.nSpoilCycles = 8;   %  number of cycles of gradient spoiling across slice thickness

nslices = 7; %20;
SLICES = [1:2:nslices 2:2:nslices];   % slice ordering (minimize slice crosstalk)

scandur = 1*60;           % seconds
delay.postrf = 10;        % (ms) delay after RF pulse. Determines TE. 
tr = 500;                 % (ms) sequence tr

fov = 22.4;          % cm
fov = 26;            % cm
nx = 64; ny = 64;    % matrix size

espmin = 0.52;      % (ms) minimum echo spacing allow by scanner (on GE, see /usr/g/bin/epiesp.dat)

gamma = 4.2576e3;   % Hz/Gauss

% number of frames
trvol = tr*nslices;  % ms
nframes = 2*ceil(scandur*1e3/trvol/2);  % force to be even

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
	'type', ex.type, ...   % 90 excitation. 'st' = small-tip
	'ftype', 'ls', ...  
	'writeModFile', false, ...
	'sliceOffset', ex.spacing, ...  % for calculating ex.freq
	'spoilDerate', 0.4, ...  % reduce slew by this much during spoiler/rewinder
	'system', system.ge);
ex.rf = toppe.utils.makeGElength(ex.rf);
ex.g = toppe.utils.makeGElength(ex.g);
toppe.writemod('rf', ex.rf, 'gx', ex.g, 'system', system.ge, 'ofname', 'tipdown.mod');

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
gyblip = toppe.utils.trapwave2(2*area/ny, system.ge.maxGrad, system.ge.maxSlew, system.ge.raster*1e3);

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

%% Calculate delay to achieve desired TR
toppe.write2loop('setup', 'version', 3);
toppe.write2loop('tipdown.mod', 'textra', delay.postrf);
toppe.write2loop('readout.mod');
toppe.write2loop('finish');
trseq = toppe.getTRtime(1,2)*1e3;    % sequence TR (ms)
delay.postreadout = max(0, tr-trseq);  % (ms) delay after readout


%% Create scanloop.txt

% rf spoiling isn't really needed since T2 << TR but why not
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;

toppe.write2loop('setup', 'version', 3);   % Initialize scanloop.txt
for ifr = 1:nframes
	% EPI calibration data:
	% frame(s)       gy    gx
	% 1, nframes-3   off   positive
	% 2, nframes-2   off   negative
	% 3, nframes-1   on    positive
	% 4, nframes     on    negative
	% turn off gy for odd/even echo calibration frames [1 2 end-1 end]
	gyamp = 1.0 - any(ifr == [1 nframes-3 2 nframes-2]);

	% flip gx for 2D odd/even echo calibration (2nd and last frame)
	gxamp = (-1)^any(ifr == [2 nframes-2 4 nframes]);

	for isl = SLICES
		% excitation
	  	toppe.write2loop('tipdown.mod', 'RFamplitude', 1.0, ...
			'RFphase', rfphs, ...
			'textra', delay.postrf, ...
			'RFoffset', round((isl-0.5-nslices/2)*ex.freq) );  % Hz (slice selection)

	 	% readout
		% data is stored in 'slice', 'echo', and 'view' indeces. Will change to ScanArchive in future.
		toppe.write2loop('readout.mod', ...
			'Gamplitude', [gxamp gyamp 0]', ...
			'DAQphase', rfphs, ...
			'slice', isl, 'echo', 1, 'view', ifr, ...  
			'textra', delay.postreadout, ...
			'dabmode', 'on');

		% update rf phase (RF spoiling)
		rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
		rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
	end
end
toppe.write2loop('finish');

tar('epi.tar', {'*.txt', '*.mod', '*.m'});

return;

% simulate and plot slice profile
fovsim = 2;          % cm
m0 = [0 0 1];        % initial magnetization
z = linspace(-fovsim/2, fovsim/2, 500);   % spatial locations (cm)
T1 = 1000;           % ms
T2 = 100;            % ms
dt = 4e-3;           % ms
figure;
[m] = toppe.utils.rf.slicesim([0 0 1], ex.rf, ex.g, dt, z, T1, T2, true);
subplot(132); title('simulated slice profile');

return;

