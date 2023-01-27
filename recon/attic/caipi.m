function IZ = caipi(ny, mb, step)
% function IZ = caipi(ny, mb, step)
% Define CAIPI sampling pattern
%
% IZ:  kz encoding index for each ky encoding location

if nargin < 3
	step = 1;
end

IZ = repmat(1:step:mb, [1 ceil(ny/mb*step)]);
IZ = IZ(1:ny);
