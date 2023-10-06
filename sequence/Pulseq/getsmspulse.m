function [rf, gz, freq] = getsmspulse(alpha, slThick, tbw, dur, nSlices, sliceSep, sysGE, sys, varargin)
% function [rf, gzRF, freq] = getsmspulse(alpha, slThick, tbw, dur, nSlices, sliceSep, sysGE, sys, varargin)
%
% Create SMS rf and gradient waveforms 
%
% Test function:
%  >> getsmspulse('test');
%
% Inputs
%   alpha        flip angle (degrees)
%   slThick      slice thickness (m)
%   tbw          time-bandwidth product
%   dur          pulse duration (sec)
%   nSlices      multi-band factor
%   sliceSep     center-to-center slice separation (m)
%   sysGE        system struct for TOPPE, see toppe.systemspecs()
%   sys          Pulseq system struct, see mr.opts()
%
% Outputs
%   rf    [nt 1]   Complex RF waveform (Gauss). Raster time is 4us.
%   g     [nt 1]   Slice-select gradient (Gauss/cm)
%   freq  [1]      Frequency offset (Hz) corresponding to sliceSep

if strcmp(alpha, 'test')
	[rf, gz] = sub_test();
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

%% Design rf and gradient (on 4us raster)
% design the 'unit' (base) pulse
nSpoilCycles = 1e-6;   % just has to be small enough so that no spoiler is added at beginning of waveform in makeslr()
[rf1, gz] = toppe.utils.rf.makeslr(alpha, slThick, tbw, dur, nSpoilCycles, sysGE, ...
	'type', arg.type, ...   
	'ftype', arg.ftype, ...  
	'ofname', arg.ofname, ...
	'writeModFile', false);

% find peak (for time reference)
I = find(abs(rf1) == max(abs(rf1(:))));
iPeak = mean(I) + 2;

% Create SMS pulse
PHS = getsmsphase(nSlices);  % Phase of the various subpulses (rad). From Wong et al, ISMRM 2012 p2209.
bw = tbw/dur*1e3;            % pulse bandwidth (Hz)
gPlateau = max(gz);           % gradient amplitude during RF excitation (Gauss/cm)
rf = 0*rf1;
dt = sysGE.raster*1e-6;      % sample (dwell) time (sec) 
t = [dt:dt:(dt*length(rf1))]' - (dt*iPeak);
for sl = 1:nSlices
	sliceOffset = (-nSlices/2 + 0.5 + sl-1) * sliceSep;   % cm
	f = sysGE.gamma*gPlateau*1e-4*sliceOffset;   % Hz
	rf = rf + rf1.*exp(1i*2*pi*f*t)*exp(1i*PHS(sl));
end

freq = sysGE.gamma*gPlateau*1e-4*sliceSep;   % Hz

%% Simulate and plot slice profile
if arg.doSim
	fov = 20;            % cm
	m0 = [0 0 1];        % initial magnetization
	z = [-fov/2:0.05:fov/2];   % spatial locations (cm)
	T1 = 1000;           % ms
	T2 = 100;            % ms
	[m] = toppe.utils.rf.slicesim([0 0 1], rf, gz, dt*1e3, z, T1, T2, true);

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


%% Create Pulseq objects

rasterIn = sysGE.raster*1e-6;    % s

% Convert from Gauss to Hz, and interpolate to sys.rfRasterTime
rfp = rf2pulseq(rf, rasterIn, sys.rfRasterTime);  

% Convert from Gauss/cm to Hz/m, and interpolate to sys.gradRasterTime
[gzp, tt] = g2pulseq(gz, rasterIn, sys.gradRasterTime);

% Trim zeros at start/end of RF pulse
% delays (s) on Siemens grad raster boundary
I = find(abs(diff(rfp) > 0));
nDelay = I(1) - mod(I(1), round(sys.gradRasterTime/sys.rfRasterTime));
rfp = rfp(nDelay:end);
delay = nDelay*sys.rfRasterTime;
I = find(abs(diff(flipud(rfp)) > 0));
rfp = rfp(1:(end-I(1)+1));

% zero-pad RF waveform duration to gradient raster boundary
wavdur = length(rfp)*sys.rfRasterTime;
ttarget = ceil(wavdur/sys.gradRasterTime)*sys.gradRasterTime;
rfp = [rfp(:); zeros(round((ttarget-wavdur)/sys.rfRasterTime), 1)];

% if delay < sys.rfDeadTime, set to rfDeadTime and delay gradients accordingly
if delay < sys.rfDeadTime
    gdelay = sys.rfDeadTime - delay;
    delay = sys.rfDeadTime;
else
    gdelay = 0;
end
    
% create pulseq objects
% Account for the fact that makeArbitraryRf scales the pulse as follows:
% signal = signal./abs(sum(signal.*opt.dwell))*flip/(2*pi);
flip = alpha/180*pi;
flipAssumed = abs(sum(rfp));
rf = mr.makeArbitraryRf(rfp, ...
    flip*abs(sum(rfp*sys.rfRasterTime))*(2*pi), ...
    'delay', delay, ...
    'system', sys);
rf.signal = rf.signal/max(abs(rf.signal))*max(abs(rfp)); % ensure correct amplitude (Hz)
gz = mr.makeArbitraryGrad('z', gzp, sys, ...
    'delay', gdelay);
gz.first = 0;
gz.last = 0;
return


function [rf, gz] = sub_test
alpha = 70;        % degrees
slThick = 5e-3;    % m
sliceSep = 40e-3;  % m
tbw = 6;
dur = 8e-3;        % s
mb = 5;            % multiband factor (number of slices)
sysGE = toppe.systemspecs();
sys = mr.opts('rfDeadTime', 100e-6, ...
    'rfRingdownTime', 60e-6, ...
    'adcRasterTime', 2e-6, ...
    'adcDeadTime', 40e-6);
[rf, gz] = getsmspulse(alpha, slThick, tbw, dur, mb, sliceSep, sysGE, sys, ...
	'doSim', true);

return;
