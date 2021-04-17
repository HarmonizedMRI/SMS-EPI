
% setup
curdir = pwd; cd ~/github/mirt; setup; cd(curdir);

% toy example
reconsms;

return;

% WIP

n = 64;  % image size

% example: MB=2, FOV/2 shift along y
ideltax = [1 1];
ideltay = [1 n/2];   % FOV/2 shift

W = getW(ideltax, ideltay, n);

L = 16;  % number of coils

N = n^2;
MB = size(ideltax, 1);  % number of simultaneous slices (multiband factor)

S_c = spalloc(N, MB*N, MB*N);
