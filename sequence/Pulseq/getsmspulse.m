function [rf, g, freq] = getsmspulse(flip, slThick, tbw, dur, nSlices, sliceSep, sys, varargin)
% function [rf, g, freq] = getsmspulse(flip, slThick, tbw, dur, nSlices, sliceSep, sys, varargin)
%
% Create SMS rf and gradient waveforms 
%
% Inputs
%   flip         flip angle (degrees)
%   slThick      slice thickness (cm)
%   tbw          time-bandwidth product
%   dur          pulse duration (msec)
%   nSlices      Multi-band factor
%   sliceSep     Center-to-center slice separation (cm)
%   sys          system struct for TOPPE, see toppe.systemspecs()
%
% Outputs
%   rf    [nt 1]   Complex RF waveform (Gauss). Raster time is 4us.
%   g     [nt 1]   Slice-select gradient (Gauss/cm)
%   freq  [1]      Frequency offset (Hz) corresponding to sliceSep

if strcmp(flip, 'test')
	sub_test();
	return;
end

slThick = slThick*100;    % cm
dur = dur*1e3;            % ms
sliceSep = sliceSep*100;  % cm

% parse inputs
arg.doSim        = false;
arg.writeModFile = false;
arg.ofname       = 'tipdown.mod';
arg.type         = 'ex';               % 'ex': 90 excitation; 'st' = small-tip
arg.ftype        = 'ls';
arg = toppe.utils.vararg_pair(arg, varargin);

% design the 'unit' (base) pulse
nSpoilCycles = 1e-3;   % just has to be small enough so that no spoiler is added at beginning of waveform in makeslr()
[rf1,g] = toppe.utils.rf.makeslr(flip, slThick, tbw, dur, nSpoilCycles, sys, ...
	'type', arg.type, ...   
	'ftype', arg.ftype, ...  
	'ofname', arg.ofname, ...
	'writeModFile', false);

% find peak (for time reference)
I = find(abs(rf1) == max(abs(rf1(:))));
iPeak = mean(I) + 2;

% Create SMS pulse
PHS = getsmsphase(nSlices);  % Phase of the various subpulses (rad). From Wong et al, ISMRM 2012 p2209.
bw = tbw/dur*1e3;          % pulse bandwidth (Hz)
gPlateau = max(g);       % gradient amplitude during RF excitation (Gauss/cm)
rf = 0*rf1;
dt = sys.raster;           % sample (dwell) time (sec) 
t = [dt:dt:(dt*length(rf1))]' - (dt*iPeak);
for sl = 1:nSlices
	sliceOffset = (-nSlices/2 + 0.5 + sl-1) * sliceSep;   % cm
	f = sys.gamma*gPlateau*sliceOffset;   % Hz
	rf = rf + rf1.*exp(1i*2*pi*f*t)*exp(1i*PHS(sl));
end

freq = sys.gamma*gPlateau*sliceSep;   % Hz

% pad to 4-sample boundary
rf = toppe.makeGElength(rf);   
g = toppe.makeGElength(g);

% simulate and plot slice profile
if arg.doSim
	fov = 20;            % cm
	m0 = [0 0 1];        % initial magnetization
	z = [-fov/2:0.05:fov/2];   % spatial locations (cm)
	T1 = 1000;           % ms
	T2 = 100;            % ms
	[m] = toppe.utils.rf.slicesim([0 0 1], rf, g, dt*1e3, z, T1, T2, true);

	% Nominal (target) slice profile (may be useful later)
	%{
	res = 0.05;          % resolution of target pattern (cm)
	n = round(fov/res);  % number of 'pixels' in target pattern
	assert( round(n/2) == n/2, 'numer of pixels in target pattern must be even');
	d = zeros(n,1);
	nSlice = slThick/res;  % number of 'pixels' across each slice
	assert(round(nSlice/2) == nSlice/2, 'number of pixels per slice in target pattern must be an even integer');
	for sl = 1:nSlices
		sliceOffset = (-nSlices/2 + 0.5 + sl-1) * sliceSep;   % cm
		nCenter = sliceOffset/res + n/2;
		assert(round(nCenter) == nCenter, "slice center location is off-grid");
		d((nCenter-nSlice/2):(nCenter+nSlice/2-1)) = 1;
	end
	subplot(1,3,2); hold on;
	x = res*(0.5:n) - res*n/2;
	plot(x,d,'g');
	legend('simulated', 'target');
	%}

end

% write to TOPPE .mod file
if arg.writeModFile
	% wrap waveforms in toppe.utils.makeGElength() to make sure they are on a 4-sample (16 us) boundary
	toppe.writemod(sys, 'rf', toppe.utils.makeGElength(rf), ...
		'gz', toppe.utils.makeGElength(g), ...
		'ofname', arg.ofname);

	% Plot .mod file. Note PNS waveform (can easily exceed 80%/100% threshold).
	%toppe.plotmod('tipdown.mod');
end

return

function sub_test
flip = 90;       % degrees
slThick = 0.5e-2;   % m
sliceSep = 4e-2;    % m
tbw = 6;
dur = 4e-3;         % s
mb = 5;          % multiband factor (number of slices)
sys = toppe.systemspecs();
getsmspulse(flip, slThick, tbw, dur, mb, sliceSep, sys, ...
	'doSim', true);
return;
