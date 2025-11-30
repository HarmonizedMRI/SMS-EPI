function [IYlabel, Jumps] = ky2blipsandrewinders(IY)
% function [IYlabel, Jumps] = ky2blipsandrewinders(IY)
%
% Categorize ky jumps as either blips (assumed small) 
% or rewinders.
% Used to determine max trapezoid area in writeEPI.m
%
% Inputs
%   IY   [etl]   ky indeces (in units of 1/fov): 1 ... ny, centered at ny/2+1
%
% Outputs
%   IYlabel   [etl-1]   true/false     true if blip, false if rewinder

etl = length(IY);

Jumps = diff(IY);   % signed jumps

% ky absolute jumps
y.absjumps = abs(diff(IY));

% we define blips as jumps < 1/2 of max jump
y.blips.indices = find(y.absjumps < max(2, max(y.absjumps)/2));

% rewinder indices (if any)
y.rewinders.indices = pge2.utils.comple(y.blips.indices, etl - 1);

% return blip/rewinder labels
IYlabel = true(1, etl-1);
IYlabel(y.rewinders.indices) = false;

