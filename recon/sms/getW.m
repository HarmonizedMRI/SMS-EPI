function getW(ideltax, ideltay, n)
%
% Image space model: I_c = (S_c*W)*x
%
% I_c    [N 1]      collapsed 2d coil image]; N = n*n; image size = [n n]
% S_c    [N MB*N]   [S_c_1' ... S_c_MB'] where S_c_l is coil sensitivity vector [N 1]
% MB     int        multiband factor (number of slices)
% W      [MB*N N]   binary mask/indexing matrix that combines aliased pixels across slices
% x      [N 1]      object to be reconstructed
% 
% ideltax   [MB 1]   x location (row index) of pixels aliasing into I_c
% ideltay   [MB 1]   y location (column index) of pixels aliasing into I_c


