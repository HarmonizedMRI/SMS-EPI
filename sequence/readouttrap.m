dt = system.ge.raster;      % sec

res = fov/nx;          % spatial resolution (cm)
kmax = 1/(2*res);      % cycles/cm
area = 2*kmax/gamma;     % G/cm * sec (area of each readout trapezoid)

% x/y prephaser
% reduce slew to reduce PNS
gpre = toppe.utils.trapwave2(area/2, system.ge.maxGrad, 0.8*system.ge.maxSlew/sqrt(2), system.ge.raster*1e3);
gpre = gpre(1:(end-1)); % remove 0 at end

% readout trapezoid
% Allow ramp sampling, and violate Nyquist slightly near kmax for now.
gxslew = 0.8*system.ge.maxSlew;  % reduce PNS
mxg = 1/(fov*gamma*dt);          % Gauss/cm
gx1 = toppe.utils.trapwave2(area, mxg, gxslew, system.ge.raster*1e3);
esp = length(gx1)*system.ge.raster*1e3;   % echo spacing (ms)
if esp > fbesp(1) & esp < fbesp(2)
	% Reduce maxGrad until echo spacing is outside forbidden range
	for s = 1:-0.02:0.1
		mxg = s*system.ge.maxGrad;
		gx1 = toppe.utils.trapwave2(area, mxg, gxslew, system.ge.raster*1e3);
		if length(gx1)*system.ge.raster*1e3 > fbesp(2)
			esp = length(gx1)*system.ge.raster*1e3; 
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
gyblip = toppe.utils.trapwave2(area/ny, system.ge.maxGrad, system.ge.maxSlew, system.ge.raster*1e3);

% gy waveform for 1st/last and other echoes
imax = find(gyblip == max(gyblip));
gyblipstart = gyblip(1:imax(1));  % first half of blip
gyblipend= gyblip(imax(2):end);
gy1 = [zeros(1,length(gx1)-length(gyblipstart)) gyblipstart]; % first echo
gyn = [gyblipend zeros(1,length(gx1)-length(gyblipend)-length(gyblipstart)) gyblipstart]; % other echoes
gylast = [gyblipend zeros(1,length(gx1)-length(gyblipend))]; % last echo

