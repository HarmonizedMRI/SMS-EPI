% Create SMS EPI fMRI scan in TOPPE (and later in Pulseq)
%
% Output:
%   TOPPE scan files:  modules.txt, scanloop.txt, tipdown.mod, and readout.mod
%
% TODO: create the corresponding Pulseq file

%% Sequence parameters

% parameters common to all sequences in this folder (SMS/2D EPI, 3D GRE coil calibration scan)
[seq, system] = getparams;
fov = seq.fov;              % in-plane field of view (cm)
nx = seq.nx; ny = seq.ny;   % matrix size
gamma = system.ge.gamma;    % Hz/Gauss
dt = system.ge.raster;      % sec

% slice-selective excitation
ex.flip = 30;        % flip angle (degrees)
ex.type = 'st';      % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
ex.ftype = 'ls';     
ex.tbw = 6;          % time-bandwidth product
ex.dur = 4;          % msec
ex.nSpoilCycles = 8;   % number of cycles of gradient spoiling across slice thickness
ex.sliceSep = seq.slThick*6;     % center-to-center separation between SMS slices (cm)
mbFactor = 6;          % sms/multiband factor (number of simultaneous slices)

nslices = 6;
if mod(nslices, mbFactor) > 0
		error('Number of slices must be multiple of MB factor');
end

% timing
delay.postrf = 1;    % (ms) delay after RF pulse. Determines TE. 
scandur = 1*20;       % seconds
tr = 300;             % (ms) If tr < minimum seq tr, minimum tr is calculated and used

SLICES = [1:2:(nslices/mbFactor) 2:2:(nslices/mbFactor)];   % slice ordering (minimize slice crosstalk)

fbesp = system.ge.forbiddenEspRange;   % ms

%% SMS excitation module
sys = system.ge;
sys.maxSlew = 8;   % G/cm/ms. Reduce PNS during slice rephaser.
[ex.rf, ex.g] = makesmspulse(ex.flip, seq.slThick, ex.tbw, ex.dur, mbFactor, ex.sliceSep, ...
	'ofname', 'tipdown.mod', ...
	'doSim', true, ...   % Plot simulated SMS slice profile
	'type', ex.type, ...
	'ftype', ex.ftype, ...
	'system', sys);

%% EPI readout module
mxs = system.ge.maxSlew;
[gx, gy, gz] = getepireadout(fov, nx, ny, ...
	system.ge.maxGrad, mxs, dt*1e3, fbesp, ...
	mbFactor, ex.sliceSep);

toppe.writemod('gx', gx, 'gy', gy, 'gz', gz, ...
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
if tr < trseq
	tr = trseq;
	fprintf('Using minimum tr (%.1f ms)\n', round(trseq));
end
	
delay.postreadout = tr-trseq;  % (ms) delay after readout

% number of frames
trvol = tr*nslices/mbFactor;  % ms
nframes = 2*ceil(scandur*1e3/trvol/2);  % force to be even

%% Create scanloop.txt
% rf spoiling isn't really needed since T2 << TR but why not
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;

toppe.write2loop('setup', 'version', 3);   % Initialize scanloop.txt
for ifr = 1:nframes
	for isl = SLICES
		% excitation
	  	toppe.write2loop('tipdown.mod', 'RFamplitude', 1.0, ...
			'RFphase', rfphs, ...
			'textra', delay.postrf); %, ...
			%'RFoffset', round((isl-0.5-nslices/2)*ex.freq) );  % Hz (slice selection)

	 	% readout
		% data is stored in 'slice', 'echo', and 'view' indeces. Will change to ScanArchive in future.
		toppe.write2loop('readout.mod', ...
			'Gamplitude', [1 1 1]', ...
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

%tar('smsepi.tar', {'*.txt', '*.mod', '*.m'});

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

