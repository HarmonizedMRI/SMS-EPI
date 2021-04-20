function IZ = caipi(ny, mb)
% Define CAIPI sampling pattern
%
% IZ = kz encoding index for each ky encoding location

IZ = repmat(1:mb, [1 ceil(ny/mb)]);
IZ = IZ(1:ny);
