% smsepi2ge.m
% Convert smsepi.seq to .tar file for execution on GE

% We set 'maxView' = etl*npartitions, such that the 'slice' index
% in the raw data file correspond to the temporal frame (image volume) index.
% The following settings must match writeSMSEPI.m
etl = 82;    % echo train length
nz = 60;     % number of slices in the reconstructed image volume
mb = 2;
npartitions = nz/mb;   % number of excitations/shots per frame

sysGE = toppe.systemspecs('maxGrad', sys.maxGrad/sys.gamma*100, ...   % G/cm
    'maxSlew', sys.maxSlew/sys.gamma/10, ...                          % G/cm/ms
    'maxRF', 0.15, ...                  % Gauss. Must be >= peak RF in sequence.
    'maxView', etl*npartitions, ...     % Determines slice/view index in data file
    'psd_rf_wait', 148, ...             % RF/gradient delay (us)
    'psd_grd_wait', 156);               % ADC/gradient delay (us)

seq2ge('smsepi.seq', sysGE, 'smsepi.tar');
