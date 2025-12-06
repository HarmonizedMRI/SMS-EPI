% convert all .seq files in current folder to .pge files

% number of ADC events for setting receive gain in auto prescan
pislquant = nz * etl;

% system hardware settings
psd_rf_wait = 100e-6;  % RF-gradient delay, scanner specific (s)
psd_grd_wait = 100e-6; % ADC-gradient delay, scanner specific (s)
b1_max = 0.25;         % Gauss
g_max = 5;             % Gauss/cm
slew_max = 20;         % Gauss/cm/ms
coil = 'xrm';          % 'hrmbuhp' (UHP); 'xrm' (MR750)
sysGE = pge2.opts(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, coil);

% pns direction weights
wt = [0.8 1 0.7];

% scan file names to convert
F = dir('*.seq');

for ii = 1:length(F)

    fn = replace(F(ii).name, {'.seq'}, '');

    if strcmp(fn, 'b0.seq'), continue, end;

    ceq = seq2ceq([fn '.seq']); 

    params = pge2.check(ceq, sysGE, 'wt', wt);

    pge2.writeceq(ceq, [fn '.pge'], 'pislquant', 10, 'params', params);

    fprintf('\n');
end

return


% plot
seq = mr.Sequence();
seq.read([fn '.seq']);
% pge2.plot(ceq, sysGE, 'blockRange', [1 10], 'rotate', false, 'interpolate', false);
pge2.validate(ceq, sysGE, seq, [], 'row', [], 'plot', false);
pge2.validate(ceq, sysGE, seq, [], 'row', [], 'plot', true);

% After simulating in WTools/VM or scanning, grab the xml files 
% and compare with the seq object:
warning('OFF', 'mr:restoreShape');  % turn off Pulseq warning for spirals
xmlPath = '~/transfer/xml/';
pge2.validate(ceq, sysGE, seq, xmlPath, 'row', [], 'plot', true);

% Check mechanical resonances (forbidden frequency bands)
S = pge2.plot(ceq, sysGE, 'blockRange', [1 10], 'rotate', true, 'interpolate', true);
check_grad_acoustics(reshape([S.gx.signal S.gy.signal S.gz.signal], [length(S.gx.signal), 1, 3])/100, 'xrm');

