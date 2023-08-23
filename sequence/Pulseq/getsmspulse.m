function [rf, gzRF, freq] = getsmspulse(alpha, slThick, tbw, dur, nSlices, sliceSep, sysGE, sys, varargin)
% function [rf, gzRF, freq] = getsmspulse(alpha, slThick, tbw, dur, nSlices, sliceSep, sysGE, sys, varargin)
%
% Create SMS rf and gradient waveforms 
%
% Inputs
%   alpha        flip angle (degrees)
%   slThick      slice thickness (m)
%   tbw          time-bandwidth product
%   dur          pulse duration (sec)
%   nSlices      Multi-band factor
%   sliceSep     Center-to-center slice separation (m)
%   sysGE        system struct for TOPPE, see toppe.systemspecs()
%   sys          Pulseq system struct, see mr.opts()
%
% Outputs
%   rf    [nt 1]   Complex RF waveform (Gauss). Raster time is 4us.
%   g     [nt 1]   Slice-select gradient (Gauss/cm)
%   freq  [1]      Frequency offset (Hz) corresponding to sliceSep

if strcmp(alpha, 'test')
	[rf, gzRF] = sub_test();
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
[rf1,g] = toppe.utils.rf.makeslr(alpha, slThick, tbw, dur, nSpoilCycles, sysGE, ...
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
gPlateau = max(g);           % gradient amplitude during RF excitation (Gauss/cm)
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


%% Create Pulseq objects

% Convert from Gauss to Hz, and interpolate to sys.rfRasterTime
wav = rf2pulseq(rf, sysGE.raster*1e-6, sys.rfRasterTime);  

% Convert from Gauss/cm to Hz/m, and interpolate to sys.gradRasterTime
gin = g;
rasterIn = sysGE.raster*1e-6;    % s
rasterOut = sys.gradRasterTime;  % s
grad.waveform = g * sysGE.gamma / 100;    % Hz/m
grad.first = grad.waveform(1);   
grad.last = grad.waveform(end);
grad.tt = (1:length(grad.waveform)) * rasterIn - rasterIn/2;
[g, tt] = gradinterp(grad, rasterIn, rasterOut);

% trim zeros at start/end of RF pulse
% delay (s) must be on grad raster boundary
I = find(abs(diff(wav) > 0));
nDelay = I(1) - mod(I(1), round(sys.gradRasterTime/sys.rfRasterTime));
wav = wav(nDelay:end);
delay = nDelay*sys.rfRasterTime;
I = find(abs(diff(fliplr(wav)) > 0));
wav = wav(1:(end-I(1)+1));

% zero-pad RF waveform duration to gradient raster boundary
wavdur = length(wav)*sys.rfRasterTime;
ttarget = ceil(wavdur/sys.gradRasterTime)*sys.gradRasterTime;
wav = [wav(:); zeros(round((ttarget-wavdur)/sys.rfRasterTime), 1)];

% if delay < sys.rfDeadTime, set to rfDeadTime and delay gradients accordingly
gdelay = sys.rfDeadTime - delay;
delay = max(sys.rfDeadTime, delay);
    
% create pulseq objects
rf = mr.makeArbitraryRf(wav, alpha/180*pi, ...
            'delay', delay, ...
            'system', sys);
gzRF = mr.makeArbitraryGrad('z', g, sys, ...
    'delay', gdelay);
gzRF.first = grad.first;

return


function [rf, gzRF] = sub_test
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
[rf, gzRF] = getsmspulse(alpha, slThick, tbw, dur, mb, sliceSep, sysGE, sys, ...
	'doSim', true);

return;
