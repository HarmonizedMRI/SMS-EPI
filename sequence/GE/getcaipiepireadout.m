function [gx, gy, gz, gpre, esp, gx1, kz, nBlipMax] = getcaipiepireadout(FOV, imSize, Ry, pf_ky, Rz, CaipiShiftZ, varargin) 
% function [gx, gy, gz, gpre, esp, gx1, kz, nBlipMax] = getcaipiepireadout(FOV, imSize, Ry, pf_ky, Rz, CaipiShiftZ, varargin) 
%
% Created 3D EPI CAIPI readout gradient waveform.
% For now, assumes isotropic resolution.
% Assumes gradient raster = ADC dwell time = 4us (for now).
%
% Inputs:
%  FOV       [1 3]  cm
%  imSize    [1 3]  image matrix size
%  Ry        [1]    ky acceleration factor
%  Rz        [1]    kz acceleration factor
%  pf_ky     [1]    Partial Fourier factor (along ky). 
%  CaipiShiftZ  [1]    size of kz step (integer multiples of 1/FOV(3))). 
% 
% Optional keyword arguments
%  gMax      [1]       nax gradient (Gauss/cm)
%  slewRead  [1 3]     max slew along x/y/z gradient axis (Gauss/cm/ms). Default: [12 15 15]
%  slewPre   [1]       max slew during prephasing gradient trapezoid (Gauss/cm/ms). Default: 12
%  fbesp     [1 2]     forbidden echo spacing range (ms). Default: [0.41 0.51]
%  plot      boolean   plot k-space? Default: false
%
% Outputs:
%  gx     x readout waveform (G/cm), without prephaser at beginning
%  gy     y readout waveform, without prephaser at beginning
%  gz     y readout waveform, without prephaser at beginning
%  gpre.x, gpre.y, gpre.z   x/y/z prephaser (to go to corner of kspace)
%  esp    echo spacing (ms)
%  gx1    waveform for one echo (G/cm)
%  kz     kz encoding indeces
%
%  See caipiepifmri.m for how to use these outputs to construct a TOPPE fMRI scan.
%
% Tips:
%  Set line 81 to 'true' to run the caipi Python call manually (just needs to create caipi.mat)

if strcmp(FOV, 'test')
    sub_test();
    return
end

if CaipiShiftZ < 1 | rem(CaipiShiftZ,1)
    error('CaipiShiftZ must be non-negative integer');
end

% Set keyword argument defaults and parse input
arg.gMax = 10;              % Peak gradient amplitude (Gauss/cm)
arg.slewRead = [11 15 15];  % Limit slew rate to this value (Gauss/cm/ms), to control PNS.
arg.slewPre = 10;   % Limit slew rate to this value (Gauss/cm/ms) during prewinder.
arg.fbesp = [0.41 0.51]; % forbidden echo spacing range (ms)
arg.plot = false;
arg.SegmentationFactor = 1;
arg.caipiPythonPath = '~/github/HarmonizedMRI/3DEPI/caipi/';

arg = toppe.utils.vararg_pair(arg, varargin);

if pf_ky < 0.7 | pf_ky > 1.0
    error('Partial Fourier factor must be in the range [0.7 1.0]');
end

nx = imSize(1); ny = imSize(2); nz = imSize(3);

dt = 4e-3;   % ms
gamma = 4257.6;        % Hz/G

res = FOV(1)/nx;       % spatial resolution (cm)
kmax = 1/(2*res);      % cycles/cm
area = 2*kmax/gamma;   % G/cm * sec (area of each readout trapezoid)
dkx = 1/FOV(1);

dky = Ry/FOV(2);       % ky spacing (cycles/cm)
dkz = 1/FOV(3);        % kz spacing (cycles/cm) corresponding to CaipiShiftZ=1

etl = ceil(pf_ky*ny/Ry);

% kz-encode indeces for one echo train
pyFile = [arg.caipiPythonPath '/skippedcaipi_sampling.py'];
pyCmd = sprintf('python %s %d %d %d %d %d %d', ...
    pyFile, ny, nz, Ry, Rz, CaipiShiftZ, 1);
fprintf('Run the following python command: %s\n', pyCmd);
if false
    input('Press any key to continue');
else
    system(pyCmd);
end
load caipi
kz1 = double(indices(:,1)' + 1);
%kz = [];
%for ii = 1:CaipiShiftZ
%    kz = [kz ii:CaipiShiftZ:Rz];
%end
%kz = repmat(kz, [1 ceil(etl/Rz)]);
kz = kz1(1:etl);

% kz encoding blip amplitudes (multiples of dkz)
kzAmp = diff(kz);

% phase-encode blips
gyBlip = toppe.utils.trapwave2(dky/gamma, arg.gMax, arg.slewRead(2), dt);  % PE1
nyBlip = length(gyBlip);
gzBlip = toppe.utils.trapwave2(dkz*max(abs(kzAmp))/gamma, arg.gMax, arg.slewRead(3), dt);  % PE2
nzBlip = length(gzBlip);
nBlipMax = max([nyBlip nzBlip]); % blip duration

% readout trapezoid (with ramp sampling)
mxg = min(1/(FOV(1)*gamma*dt*1e-3), arg.gMax);    % Gauss/cm
gx1 = toppe.utils.trapwave2(area, mxg, arg.slewRead(1), dt);

% Extend readout plateau to make room for turns
gx1orig = gx1;
nExtend = 0; % number of samples to add in middle of gx1
area1 = sum(gx1(ceil(nBlipMax/2):(end-ceil(nBlipMax/2))))*dt*1e-3;   % G/cm/s
while area1 < area
    nExtend = nExtend + 1;
    gx1 = [gx1orig(1:ceil(end/2)) max(gx1)*ones(1,nExtend) gx1orig((ceil(end/2)+1):end)];
    area1 = sum(gx1(ceil(nBlipMax/2):(end-ceil(nBlipMax/2))))*dt*1e-3;   % G/cm/s
end

% Reduce peak gradient until echo spacing is outside forbidden range
esp = length(gx1)*dt;  % echo spacing (ms)
if esp > arg.fbesp(1) & esp < arg.fbesp(2)
    for s = 1:-0.005:0.1
        mxg = s*mxg;
        gx1 = toppe.utils.trapwave2(area, mxg, arg.slewRead(1), dt);
        if length(gx1)*dt > arg.fbesp(2)
            esp = length(gx1)*dt;
            break;
        end
    end
end

% remove 0 at end in preparation for assembling into echo train
gx1 = gx1(1:(end-1));  

% echo train
gx = [];
for iecho = 1:etl
    gx = [gx gx1*(-1)^(iecho+1)];
end
gx = [gx(:); 0];  % waveform must end with 0

% add y and z blips
nt = length(gx1);  % number of samples in one echo
gy = 0*gx;
gz = 0*gx;
for ii = 1:(etl-1)
    iStart = ii*nt - nyBlip/2;
    iStop = iStart + nyBlip-1;
    gy(iStart:iStop) = gyBlip;

    iStart = ii*nt - nzBlip/2;
    iStop = iStart + nzBlip-1;
    gz(iStart:iStop) = gzBlip * kzAmp(ii)/max(abs(kzAmp));
end

% make length multiple of 4 
gx = toppe.makeGElength(gx(:));
gy = toppe.makeGElength(gy(:));
gz = toppe.makeGElength(gz(:));

% x/y/z prephasers.
areax = sum(gx1)*dt*1e-3;   % G/cm * sec
gpre.x = toppe.utils.trapwave2(areax/2, arg.gMax, arg.slewPre, dt);
gpre.y = toppe.utils.trapwave2(area/2-dky/Ry/gamma/2, arg.gMax, arg.slewPre, dt);
gpre.y = [gpre.y(:); zeros(length(gpre.x)-length(gpre.y), 1)]; % make same length
gpre.x = toppe.makeGElength(gpre.x(:)); % make length multiple of 4
gpre.y = toppe.makeGElength(gpre.y(:));
gpre.z = gpre.y; % isotropic resolution

if arg.plot
    % plot k-space. Add prephaser for plotting purposes only.
    figure;
    gxf = [-gpre.x; gx];
    gyf = [-gpre.y; gy];
    gzf = [-gpre.z; gz];
    kxp = gamma*dt*1e-3*cumsum(gxf);
    kyp = gamma*dt*1e-3*cumsum(gyf);
    kzp = gamma*dt*1e-3*cumsum(gzf);
    subplot(121); plot(kxp,kyp,'b.'); axis equal;
    xlabel('kx (cycles/cm)'); ylabel('ky (cycles/cm)');
    T = dt*(1:length(kxp));
    subplot(122); hold off; plot(T,kxp,'r'); hold on; plot(T,kyp,'g'); plot(T,kzp,'b'); hold off;
    legend('kx', 'ky', 'kz'); ylabel('cycles/cm'); xlabel('time (ms)');
end

return


function sub_test()

    res = 0.275;  % cm
    imSize = [80 80 60];
    FOV = imSize*res;  % cm
    Ry = 1;
    pf_ky = 0.8;
    Rz = 6;
    CaipiShiftZ = 2;

    [gx, gy, gz, gpre, esp, gx1, kz, nBlipMax] = getcaipiepireadout(FOV, imSize, Ry, pf_ky, Rz, CaipiShiftZ, ...
        'plot', true) ;

    load caipi
    figure; im(mask);  % ky/kz sampling pattern

return
