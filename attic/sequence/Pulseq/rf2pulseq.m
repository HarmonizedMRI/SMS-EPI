function [rfOut, ttOut] = rf2pulseq(rf, rasterIn, rasterOut)
% function [rfOut, ttOut] = rf2pulseq(rf, rasterIn, rasterOut)
%
% Convert rf waveform from Gauss to Hz, and interpolate to rasterOut
%
% Inputs
%   rf            [n 1]   RF waveform, Gauss
%   rasterIn      [1]     Input waveform raster time (sec)
%   rasterOut     [1]     Output waveform raster time (sec)
%
% Output
%   rfOut         [m 1]   RF waveform, Hz
%   ttOut         [m 1]   Output waveform sample times (sec)

gamma = 4.2576e3;    % Hz/G
rf = rf*gamma;       % Hz

dur = numel(rf)*rasterIn;   % pulse duration
ttIn = (1:length(rf))*rasterIn - rasterIn/2;
ttOut = rasterOut/2:rasterOut:dur;

rfOut = interp1(ttIn, rf, ttOut, 'linear', 'extrap');

% return column vectors
rfOut = rfOut(:);   
ttOut = ttOut(:);   

%L = 10; cutoff = 0.9;
%rf = interp(rf,dt/rfRasterTime,L,cutoff);      % upsample from 4us to 1us
