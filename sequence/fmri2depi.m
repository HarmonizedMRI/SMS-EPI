% Create 2D EPI fMRI scan in TOPPE (and later in Pulseq)
%
% Output:
%   TOPPE scan files:  modules.txt, scanloop.txt, tipdown.mod, and readout.mod
%
% TODO: create the corresponding Pulseq file

isCalScan = false;    % measure kspace and B0 eddy current using Duyn's method 

%% Sequence parameters

% parameters common to all sequences in this folder (SMS/2D EPI, 3D GRE coil calibration scan)
[seq, system] = getparams;
fov = seq.fov;              % cm
nx = seq.nx; ny = seq.ny;   % matrix size
gamma = system.ge.gamma;    % Hz/Gauss

% slice-selective excitation
ex.flip = 45;        % flip angle (degrees)
ex.type = 'st';      % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
ex.ftype = 'ls';     
ex.tbw = 6;          % time-bandwidth product
ex.dur = 4;          % msec
ex.nSpoilCycles = 8;   %  number of cycles of gradient spoiling across slice thickness

delay.postrf = 1;        % (ms) delay after RF pulse. Determines TE. 

if isCalScan
	ex.thick = 0.3;    % slice thickness (cm)
	ex.spacing = 10;   % center-to-center slice separation (cm)
	nslices = 7;
	scandur = 30;      % seconds
	tr = 500;          % (ms) sequence tr
else
	ex.thick = seq.slThick;    % slice thickness (cm)
	ex.spacing = ex.thick;  % center-to-center slice separation (cm)
	nslices = nx;      
	scandur = 1*60;    % seconds
	tr = 100;            % (ms) If tr < minimum seq tr, minimum tr is calculated and used
end

SLICES = [1:4:nslices 2:4:nslices 3:4:nslices 4:4:nslices];   % slice ordering (minimize slice crosstalk)

fbesp = system.ge.forbiddenEspRange;   % ms

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
if isCalScan 
	% excite slice along x (readout direction) so we can measure eddy currents
	toppe.writemod('rf', ex.rf, 'gx', ex.g, 'system', system.ge, 'ofname', 'tipdown.mod');
else
	toppe.writemod('rf', ex.rf, 'gz', ex.g, 'system', system.ge, 'ofname', 'tipdown.mod');
end

%% EPI readout
readouttrap;  % creates gpre and gx1. Also used in fmri2depi.m

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
if tr < trseq
	tr = trseq;
	fprintf('Using minimum tr (%.1f ms)\n', round(trseq));
end
	
delay.postreadout = tr-trseq;  % (ms) delay after readout

% number of frames
trvol = tr*nslices;  % ms
nframes = 2*ceil(scandur*1e3/trvol/2);  % force to be even


%% Create scanloop.txt

% rf spoiling isn't really needed since T2 << TR but why not
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;

toppe.write2loop('setup', 'version', 3);   % Initialize scanloop.txt
for ifr = 1:nframes
	% Calibration data:
	% frame(s)       gy    gx
	% 1, nframes-4   off   positive
	% 2, nframes-3   off   negative
	% 3, nframes-2   on    positive
	% 4, nframes-1   on    negative
	% 5, nframes     off   off

	% turn off gy for odd/even echo calibration frames
	gyamp = 1.0 - any(ifr == [1 nframes-4 2 nframes-3 5 nframes]);

	% flip gx for odd/even echo calibration
	% turn off for frames [5 nframes]
	gxamp = (-1)^any(ifr == [2 nframes-3 4 nframes-1]);
	gxamp = gxamp*(1.0 - any(ifr == [5 nframes]));

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

