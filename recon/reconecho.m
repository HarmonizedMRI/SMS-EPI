function [x,A,dcf] = reconecho(y, A, dcf, kx, nx, fov, gx)
% 1D nufft reconstruction
%
% y:   [nt]   
% kx:  [nt]  (cycles/cm)
% nx   int
% fov  cm
% gx   [nt]  (a.u.) readout gradient (typically a trapezoid) for one echo in EPI train

if strcmp(y, "test")
	sub_test;
	return
end


if isempty(A)
	nufft_args = {[nx],[6],[2*nx],[nx/2],'minmax:kb'};
	mask = true(nx,1);
	A = Gmri([fov*kx(:)],mask,'nufft',nufft_args);
end
if isempty(dcf)
	dcf = gx(:)/max(gx)/nx;  % density compensation
end

x = A'*(y(:).*dcf(:));

return


function sub_test

% 1d test object
nx = 128;
p = phantom(nx);
x = p(:,end/2);
fov = 20;  % cm

% nonuniform sampling (trapezoidal gradient)
dt = 4e-6;             % gradient raster time (s)
gamma = 4257.6;        % Hz/G
res = fov/nx;          % spatial resolution (cm)
kmax = 1/(2*res);      % cycles/cm
area = kmax/gamma;     % G/cm * sec (area of each readout trapezoid)
gmax = 1/(fov*gamma*dt);    % Gauss/cm
gslew = 10;      % G/cm/ms
gx = toppe.utils.trapwave2(2*area, gmax, gslew, dt*1e3);
gx = gx(2:(end-1));
kx = gamma*dt*cumsum(gx);
kx = kx - max(kx)/2;

% data
nufft_args = {[nx],[6],[2*nx],[nx/2],'minmax:kb'};
mask = true(nx,1);
A = Gmri([fov*kx(:)],mask,'nufft',nufft_args);
y = A*x(:);

% recon
[~,A,dcf] = reconecho(y, [], [], kx, nx, fov, gx);
xhat = reconecho(y, A, dcf); %, kx, nx, fov);
hold off; plot(x); hold on; plot(abs(xhat),'o'); legend('x true', 'xhat');

return;


kmax = max(kx);
datg = interp1(kx, dat, linspace(kx(1),kx(end),nx));
if kx(1) > kx(end)
	datg = flipdim(datg,2);
end
datg(isnan(datg)) = 0;
x = fftshift(ifft(fftshift(datg(:),1), [], 1),1);
return;
