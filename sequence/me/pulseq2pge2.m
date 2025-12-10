function pulseq2pge2(fn)
% convert .seq file to .pge file
% input:
%   fn: seq-file name

% convert to .pge
ceq = seq2ceq([fn '.seq']);   %, 'usesRotationEvents', false);

psd_rf_wait = 100e-6;  % RF-gradient delay, scanner specific (s)
psd_grd_wait = 100e-6; % ADC-gradient delay, scanner specific (s)
b1_max = 0.25;         % Gauss
g_max = 5;             % Gauss/cm
slew_max = 20;         % Gauss/cm/ms
coil = 'xrm';          % 'hrmbuhp' (UHP); 'xrm' (MR750)
sysGE = pge2.opts(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, coil);

wt = 1*[0.8 1 0.7];
%wt = [1 1 1];
params = pge2.check(ceq, sysGE, 'wt', wt);

pge2.writeceq(ceq, [fn '.pge'], 'pislquant', 1, 'params', params);

% plot
seq = mr.Sequence();
seq.read([fn '.seq']);
pge2.plot(ceq, sysGE, 'blockRange', [1 10], 'rotate', false, 'interpolate', false);
%pge2.validate(ceq, sysGE, seq, [], 'row', [], 'plot', false);
pge2.validate(ceq, sysGE, seq, [], 'row', [], 'plot', true);
end