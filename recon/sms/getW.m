function W = getW(ideltax, ideltay, n)
%
% n         [1 1]    int. image matrix size is [n n]
% ideltax   [MB 1]   x location (row index) of pixels aliasing into I_c
% ideltay   [MB 1]   y location (column index) of pixels aliasing into I_c

MB = size(ideltax, 1);  % number of simultaneous slices (multiband factor)

N = n^2;

W = spalloc(MB*N, N, N);
IJ = sub2ind([n n], ideltax, ideltay)
W(IJ(1)
for isl = 1:MB
	i = IJ(isl)
	W(

