function IZ = caipi(ny, mb, skip)
% Define CAIPI sampling pattern
%
% IZ = kz encoding index for each ky encoding location

IZ = repmat(1:skip:mb, [1 ceil(ny/mb*skip)]);
IZ = IZ(1:ny);
