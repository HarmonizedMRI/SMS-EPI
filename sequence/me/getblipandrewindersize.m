function [maxBlip, maxRewinder] = getblipandrewindersize(IY)
% function [maxBlip, maxRewinder] = getblipandrewindersize(IY)
%
% Categorize ky jumps as either 'blips' (assumed small) 
% or 'rewinders', and get corresponding max k-space jumps.
% Used to determine max trapezoid area in writeEPI.m
%
% Inputs
%   IY   [etl]   ky indeces (in units of 1/fov): 1 ... ny, centered at ny/2+1
%
% Outputs
%   maxBlip       [1]    max blip size in ky [1/fov]
%   maxRewinder   [1]    max rewinder size in ky [1/fov]

% ky jumps
y.jumps = abs(diff(IY));

% we defin blips as jumps < 1/2 of max jump
y.blips.indices = find(y.jumps < max(1, max(y.jumps)/2));

% max blip size
maxBlip = max(y.jumps(y.blips.indices));

% max rewinder size 
y.rewinders.indices = comple(y.blips.indices, length(y.jumps));
maxRewinder = pge2.utils.iff(isempty(y.rewinders.indices), [], max(y.jumps(y.rewinders.indices)));

