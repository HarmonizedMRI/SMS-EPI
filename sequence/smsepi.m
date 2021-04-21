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
ex.flip = 45;        % flip angle (degrees)
ex.type = 'st';      % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
ex.ftype = 'ls';     
ex.tbw = 6;          % time-bandwidth product
ex.dur = 4;          % msec
ex.nSpoilCycles = 8;   % number of cycles of gradient spoiling across slice thickness
ex.sliceSep = 5;       % center-to-center separation between SMS slices (cm)
mbFactor = 4;          % sms/multiband factor (number of simultaneous slices)

nslices = 3;
if mod(nslices, mbFactor) > 0
		error('Number of slices must be multiple of MB factor');
end

% timing
delay.postrf = 1;    % (ms) delay after RF pulse. Determines TE. 
scandur = 1*20;       % seconds
tr = 500;             % (ms) If tr < minimum seq tr, minimum tr is calculated and used

SLICES = [1:2:(nslices/mbFactor) 2:2:(nslices/mbFactor)];   % slice ordering (minimize slice crosstalk)

fbesp = system.ge.forbiddenEspRange;   % ms

%% SMS excitation module
sys = system.ge;
sys.maxSlew = 8;   % G/cm/ms. Reduce PNS during slice rephaser.
if 0
[ex.rf, ex.g] = makesmspulse(ex.flip, seq.slThick, ex.tbw, ex.dur, mbFactor, ex.sliceSep, ...
	'ofname', 'tipdown.mod', ...
	'doSim', true, ...   % Plot simulated SMS slice profile
	'system', sys);
end

%% EPI readout
res = fov/nx;          % spatial resolution (cm)
kmax = 1/(2*res);      % cycles/cm
area = 2*kmax/gamma;     % G/cm * sec (area of each readout trapezoid)

% x/y prephaser
% reduce slew to reduce PNS
gpre = toppe.utils.trapwave2(area/2, system.ge.maxGrad, 0.8*system.ge.maxSlew/sqrt(2), system.ge.raster*1e3);
gpre = gpre(1:(end-1)); % remove 0 at end

% readout trapezoid
% Allow ramp sampling, and violate Nyquist slightly near kmax for now.
gxslew = 0.8*system.ge.maxSlew;  % reduce PNS
mxg = 1/(fov*gamma*dt);          % Gauss/cm
gx1 = toppe.utils.trapwave2(area, mxg, gxslew, system.ge.raster*1e3);
esp = length(gx1)*system.ge.raster*1e3;   % echo spacing (ms)
if esp > fbesp(1) & esp < fbesp(2)
	% Reduce maxGrad until echo spacing is outside forbidden range
	for s = 1:-0.02:0.1
		mxg = s*system.ge.maxGrad;
		gx1 = toppe.utils.trapwave2(area, mxg, gxslew, system.ge.raster*1e3);
		if length(gx1)*system.ge.raster*1e3 > fbesp(2)
			esp = length(gx1)*system.ge.raster*1e3; 
			break;
		end
	end
end
gx1 = gx1(1:(end-1));  % remove 0 at end

% check that Nyquist is supported everywhere along readout
kx = gamma*dt*cumsum(gx1);
minfov = 1/max(abs(diff(kx)));
if minfov < fov
	error('Nyquist violated along readout direction');
end

% y blip
gyblip = toppe.utils.trapwave2(area/ny, system.ge.maxGrad, system.ge.maxSlew, system.ge.raster*1e3);

% gy waveform for 1st/last and other echoes
imax = find(gyblip == max(gyblip));
gyblipstart = gyblip(1:imax(1));  % first half of blip
gyblipend= gyblip(imax(2):end);
gy1 = [zeros(1,length(gx1)-length(gyblipstart)) gyblipstart]; % first echo
gyn = [gyblipend zeros(1,length(gx1)-length(gyblipend)-length(gyblipstart)) gyblipstart]; % other echoes
gylast = [gyblipend zeros(1,length(gx1)-length(gyblipend))]; % last echo

% z blip/rewinder. Use one waveform for both and scale as needed.
kmax = 1/(2*ex.sliceSep);   % cycles/cm
area = 2*kmax/gamma;        % G/cm * sec
gzblip = toppe.utils.trapwave2(area, system.ge.maxGrad, system.ge.maxSlew, system.ge.raster*1e3);

% z prewinder
gzpre = [(mbFactor/2-1/2)/2*gzblip zeros(1,length(gpre)-length(gzblip))];

% gz waveforms for the various echoes
imax = find(gzblip == max(gzblip));
gzblipstart = gzblip(1:imax(1));  % first half of blip
gzblipend = gzblip(imax(2):end);
amp = 1/(mbFactor-1);  % amplitude of one delta_kz step (scale gzblip by amp)
gz1 = [zeros(1,length(gx1)-length(gzblipstart)) amp*gzblipstart]; % first echo
gz2 = [amp*gzblipend zeros(1,length(gx1)-length(gzblipend)-length(gzblipstart)) amp*gzblipstart];
gz3 = [amp*gzblipend zeros(1,length(gx1)-length(gzblipend)-length(gzblipstart)) -gzblipstart];
gz4 = [-gzblipend zeros(1,length(gx1)-length(gzblipend)-length(gzblipstart)) amp*gzblipstart];
gzlast = [amp*gzblipend zeros(1,length(gx1)-length(gzblipend))]; 

% put it all together 
gx = [-gpre gx1];
gy = [-gpre gy1];
gz = [-gzpre gz1];
for iecho = 2:(ny-1)
	gx = [gx gx1*(-1)^(iecho+1)];
	gy = [gy gyn];
end
gx = [gx gx1*(-1)^(iecho+2) 0];  % add zero at end
gy = [gy gylast 0];
for iecho=2:mbFactor:(ny-mbFactor+1)
	gz = [gz repmat(gz2, 1, mbFactor-2) gz3 gz4];
end
gz = [gz gzlast 0];

% fix gz at end
gz = gz(1:length(gx));
gz((end-round(length(gzlast)/2)+1):end) = 0;

% write to readout.mod
gx = toppe.utils.makeGElength(gx(:));
gy = toppe.utils.makeGElength(gy(:));
gz = toppe.utils.makeGElength(gz(:));
toppe.writemod('gx', gx, 'gy', gy, 'gz', gz, ...
	'system', system.ge, ...
	'ofname', 'readout.mod');

% plot kspace
[~,gx,gy,gz] = toppe.readmod('readout.mod');
kx = gamma*dt*cumsum(gx);
ky = gamma*dt*cumsum(gy);
kz = gamma*dt*cumsum(gz);
figure;
subplot(121); plot(kx,ky,'b.'); axis equal;
xlabel('kx (cycles/cm)'); ylabel('ky (cycles/cm)');
T = dt*1e3*(1:length(kx));
subplot(122); hold off; plot(T,kx,'r'); hold on; plot(T,ky,'g'); plot(T,kz,'b'); hold off;
legend('kx', 'ky', 'kz'); ylabel('cycles/cm'); xlabel('time (ms)');


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

tar('smsepi.tar', {'*.txt', '*.mod', '*.m'});

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

