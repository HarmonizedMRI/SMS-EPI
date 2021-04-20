function [x,A,dcf] = reconecho(y, nx, A, dcf, kx, fov)
% function [x,A,dcf] = reconecho(y, nx, [A, dcf], [kx, fov])
%
% 1D nufft reconstruction of ramp-sampled EPI echo
%
% y:    [nt]  If empty, x is empty.
% nx    int
% A:    Gmri object. If empty, construct from kx, nx, fov
% dcf:  [nt] density compensation. If empty, construct from kx.
% kx:   [nt]  (cycles/cm)
% fov   cm

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
	gx = abs([diff(kx(:)); 0]);
	dcf = gx(:)/max(gx);
end

if ~isempty(y)
	x = A'*(y(:).*dcf(:))/nx;
else
	x = [];
end

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
% get A and dcf first, then apply (way faster)
[~,A,dcf] = reconecho([], nx, [], [], kx, fov);
nops = 500;
fprintf('Doing %d 1d recons... ', nops);
tic;
for ii = 1:nops
	xhat = reconecho(y, nx, A, dcf);
end
toc;
hold off; plot(x); hold on; plot(abs(xhat),'o'); legend('x true', 'xhat');

if 0
% compare with fft
kx = linspace(kx(1), kx(end), nx);
[~,A,dcf] = reconecho([], nx, [], [], kx, fov);
y = A*x;
yf = fftshift(fft(fftshift(x)));
figure; plot(abs(y)); hold on; plot(abs(yf)); legend('y', 'ift');

x1 = A'*y;
xf = fftshift(ifft(fftshift(yf)));
figure; hold on; plot(x); plot(abs(x1)); plot(abs(xf),'o'); legend('x true', 'xhat', 'fft');
end

return;


kmax = max(kx);
datg = interp1(kx, dat, linspace(kx(1),kx(end),nx));
if kx(1) > kx(end)
	datg = flipdim(datg,2);
end
datg(isnan(datg)) = 0;
x = fftshift(ifft(fftshift(datg(:),1), [], 1),1);
return;
