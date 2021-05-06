function [gx, gy, gz] = getepireadout(fov, nx, ny, mb, sliceSep, gMax, slewMax, raster, fbesp)

dt = raster*1e-3;      % sec
gamma = 4257.6;        % Hz/G

res = fov/nx;          % spatial resolution (cm)
kmax = 1/(2*res);      % cycles/cm
area = 2*kmax/gamma;     % G/cm * sec (area of each readout trapezoid)

% x/y prephaser
% reduce slew to reduce PNS
gpre = toppe.utils.trapwave2(area/2, gMax, 0.8*slewMax/sqrt(2), dt*1e3);
gpre = gpre(1:(end-1)); % remove 0 at end

% readout trapezoid
% Allow ramp sampling, and violate Nyquist slightly near kmax for now.
gxslew = 0.8*slewMax;  % reduce PNS
mxg = 1/(fov*gamma*dt);          % Gauss/cm
gx1 = toppe.utils.trapwave2(area, mxg, gxslew, dt*1e3);
esp = length(gx1)*dt*1e3;   % echo spacing (ms)
if esp > fbesp(1) & esp < fbesp(2)
	% Reduce maxGrad until echo spacing is outside forbidden range
	for s = 1:-0.02:0.1
		mxg = s*gMax;
		gx1 = toppe.utils.trapwave2(area, mxg, gxslew, dt*1e3);
		if length(gx1)*dt*1e3 > fbesp(2)
			esp = length(gx1)*dt*1e3; 
			break;
		end
	end
end
gx1 = gx1(1:(end-1));  % remove 0 at end

% check that Nyquist is supported everywhere along readout
kx = gamma*dt*cumsum(gx1);
minfov = 1/max(abs(diff(kx)));
if minfov < fov
	error('Nyquist violated along readout direction');
end

% y blip
gyblip = toppe.utils.trapwave2(area/ny, gMax, slewMax, dt*1e3);

% gy waveform for 1st/last and other echoes
imax = find(gyblip == max(gyblip));
gyblipstart = gyblip(1:imax(1));  % first half of blip
gyblipend= gyblip(imax(2):end);
gy1 = [zeros(1,length(gx1)-length(gyblipstart)) gyblipstart]; % first echo
gyn = [gyblipend zeros(1,length(gx1)-length(gyblipend)-length(gyblipstart)) gyblipstart]; % other echoes
gylast = [gyblipend zeros(1,length(gx1)-length(gyblipend))]; % last echo

% z blip/rewinder. Use one waveform for both and scale as needed.
kmax = 1/(2*sliceSep);   % cycles/cm
area = 2*kmax/gamma;        % G/cm * sec
gzblip = toppe.utils.trapwave2(area, gMax, slewMax, dt*1e3);

% z prewinder
gzpre = [(mb/2-1/2)/2*gzblip zeros(1,length(gpre)-length(gzblip))];
gzpre = [gzblip/2 zeros(1,length(gpre)-length(gzblip))];

% gz waveforms for the various echoes
imax = find(gzblip == max(gzblip));
gzblipstart = gzblip(1:imax(1));  % first half of blip
gzblipend = gzblip(imax(2):end);
amp = 1/(mb/2-1);  % amplitude of one delta_kz step (scale gzblip by amp)
gz1 = [zeros(1,length(gx1)-length(gzblipstart)) amp*gzblipstart]; % first echo
gz2 = [amp*gzblipend zeros(1,length(gx1)-length(gzblipend)-length(gzblipstart)) amp*gzblipstart];
gz3 = [amp*gzblipend zeros(1,length(gx1)-length(gzblipend)-length(gzblipstart)) -gzblipstart];
gz4 = [-gzblipend zeros(1,length(gx1)-length(gzblipend)-length(gzblipstart)) amp*gzblipstart];
gzlast = [amp*gzblipend zeros(1,length(gx1)-length(gzblipend))]; 

% put it all together 
% include 3 calibration echoes at start
gx = [gpre -gx1 gx1 -gx1];
gy = 0*gx;
gz = 0*gx;
gx = [gx 0*gpre  gx1];
gy = [gy  -gpre  gy1];
gz = [gz  -gzpre gz1];
for iecho = 2:(ny-1)
	gx = [gx gx1*(-1)^(iecho+1)];
	gy = [gy gyn];
end
gx = [gx gx1*(-1)^(iecho+2) 0];  % add zero at end
gy = [gy gylast 0];
for iecho = 2:mb/2:ny
	gz = [gz repmat(gz2, 1, mb/2-2) gz3 gz4];
end
gz = [gz gzlast 0];

% fix gz at end
gz = gz(1:length(gx));
gz((end-round(length(gzlast)/2)+1):end) = 0;

gx = toppe.utils.makeGElength(gx(:));
gy = toppe.utils.makeGElength(gy(:));
gz = toppe.utils.makeGElength(gz(:));

% plot kspace
kx = gamma*dt*cumsum(gx);
ky = gamma*dt*cumsum(gy);
kz = gamma*dt*cumsum(gz);
figure;
subplot(121); plot(kx,ky,'b.'); axis equal;
xlabel('kx (cycles/cm)'); ylabel('ky (cycles/cm)');
T = dt*1e3*(1:length(kx));
subplot(122); hold off; plot(T,kx,'r'); hold on; plot(T,ky,'g'); plot(T,kz,'b'); hold off;
legend('kx', 'ky', 'kz'); ylabel('cycles/cm'); xlabel('time (ms)');

