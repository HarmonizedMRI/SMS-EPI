function [gOut, ttOut] = g2pulseq(g, rasterIn, rasterOut)
% function [gOut, ttOut] = g2pulseq(g, rasterIn, rasterOut)
%
% Convert grad waveform from G/cm to Hz/m, and interpolate to rasterOut
%
% Inputs
%   g            [n 1]   gradient waveform, Gauss/cm
%   rasterIn     [1]     Input waveform raster time (sec)
%   rasterOut    [1]     Output waveform raster time (sec)
%
% Output
%   gOut         [m 1]   gradient waveform, Hz/m
%   ttOut        [m 1]   Output waveform sample times (sec)

gamma = 4.2576e3;    % Hz/G
g = g * gamma * 100;   % Hz/m

if rasterOut < rasterIn
    error('rasterOut must be > rasterIn');
end

ttIn = (1:length(g))*rasterIn - rasterIn/2;
ttOut = rasterOut/2:rasterOut:ttIn(end);

gOut = interp1(ttIn, g, ttOut); 

if any(isnan(gOut))
    error('g2pulseq(): NaN in gradient waveform after interpolation');
end

% return column vectors
gOut = gOut(:);   
ttOut = ttOut(:);   
