
% Get odd/even echo correction parameters 
fn = [datdir pfileCal]; 
ifr = 1;
draw = hmriutils.epi.loadframeraw_ge(fn, etl, np, ifr);   % [nFID etl np nc]

%[d, mask] = readpartition(fn, N, mb, etl, frame, round(np/2), sysGE, fovXcm, [0 0]', readoutFile, caipiFile);
nFID = size(d, 1);                             % number of data samples in ADC window

% get calibration factor a
ds = sum(d, 3);
a = getcalparams(squeeze(ds));
save a a

% Get single-slice EPI images
main1;


