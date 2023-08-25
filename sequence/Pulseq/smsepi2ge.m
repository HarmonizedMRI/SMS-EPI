% smsepi2ge.m
% Convert smsepi.seq to .tar file for execution on GE

sysGE = toppe.systemspecs('maxGrad', sys.maxGrad/sys.gamma*100, ...   % G/cm
    'maxSlew', sys.maxSlew/sys.gamma/10, ...                          % G/cm/ms
    'maxRF', 0.15, ...                  % Gauss. Must be >= peak RF in sequence.
    'maxView', 500, ...              % Determines slice/view index in data file
    'psd_rf_wait', 148, ...          % RF/gradient delay (us)
    'psd_grd_wait', 156);            % ADC/gradient delay (us)

seq2ge('smsepi.seq', sysGE, 'smsepi.tar');
