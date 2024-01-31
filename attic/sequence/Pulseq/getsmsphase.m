function PHS = getSMSphase(mb)
% Return phase of SMS subpulses (up to mb=8) from Table 1 in Wong ISMRM 2012 p2209
%
% Input
%   mb    multiband factor (1-8, int)
%
% Output
%   PHS   [1 mb]  radians

mbMax = 8;
if mb > mbMax
	error(sprintf('Max SMS factor is %d', mbMax));
end

PHS = zeros(8,8);
PHS(1,1) = 0;
PHS(2,1:2) = [0 pi];  % arbitrary
PHS(3,1:3) = [0 0.730 4.602];
PHS(4,1:4) = [0 3.875 5.940 6.197];
PHS(5,1:5) = [0 3.778 5.335 0.872 0.471];
PHS(6,1:6) = [0 2.005 1.674 5.012 5.736 4.123];
PHS(7,1:7) = [0 3.002 5.998 5.909 2.624 2.528 2.440];
PHS(8,1:8) = [0 1.036 3.414 3.778 3.215 1.756 4.555 2.467];

PHS = PHS(mb,1:mb);

return;
