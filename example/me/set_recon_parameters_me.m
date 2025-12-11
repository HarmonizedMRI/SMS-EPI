% set recon parameters for multi-echo sms epi
nFrames = 40;
nFramesDiscard = 12;
%ncoil = 32;

iecho =2; % which echo to be recon'd

ikz = [1 3 1 3];
Ry = 2; % in-plane acce factor

% calibration region for slice-grappa
ncalx = 48;
ncaly = 44;
Calx = nx/2-ncalx/2:nx/2+ncalx/2-1;
Caly = ny/2-ncaly/2:ny/2+ncaly/2-1;
Caly = Caly - (ny-etl*Ry);


% calibration region for in-plane grappa
ncalx = 48;
ncaly = 48;
Calx_grappa = nx/2 - ncalx/2:nx/2+ncalx/2-1;
Caly_grappa = ny/2 - ncaly/2:ny/2+ncaly/2-1;
