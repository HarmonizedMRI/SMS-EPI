function [IY, IZ] = getcaipi(ny, nz, Ry, Rz, caipiShiftZ, pth)
% function [IY, IZ] = getcaipi(ny, nz, Ry, Rz, caipiShiftZ, skippedCaipiPath)
%
% Get CAIPI sampling pattern (for one shot/echo train)
%
% Output
%   IY   [etl]   ky indeces (in units of 1/fov): 1 ... ny, centered at ny/2+1
%   IZ   [etl]   kz indeces (in units of 1/fov): 1 ... Rz, centered at Rz/2+1

pth = pge2.utils.normalizepath(pth);

assert(~mod(ny,2), 'ny must be an even integer');
assert(~mod(nz,2), 'ny must be an even integer');
assert(~mod(Ry,1), 'Ry must be an integer');
assert(~mod(Rz,1), 'Rz must be an integer');
assert(~mod(caipiShiftZ,1), 'caipiShiftZ must be an integer');

% create caipi.mat, then load it
pyFile = [pth 'skippedcaipi_sampling.py'];
pyCmd = sprintf('python3 %s %d %d %d %d %d %d', ...
    pyFile, ny, nz, Ry, Rz, caipiShiftZ, 1);

if system(pyCmd) ~= 0
    fprintf('Open a terminal and run the following python command:\n\t%s\n', pyCmd);
    input('\tWhen done, press Enter to continue');
end

load caipi

% kz and ky indeces (multiples of deltak)
IY = double(indices(:, 2)) + 1;
IZ = double(indices(:, 1)) + 1;
