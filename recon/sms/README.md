
Image space model: 

```
I_c = (S_c*W)*x    [N 1]      collapsed 2d coil image]; N = n*n; image size = [n n]

S_c    [N MB*N]   S_c(1,:) = [S_c_1' ... S_c_MB'] where S_c_l is coil sensitivity map [N 1]
                  Each row is identical.
W      [MB*N N]   indicator/indexing matrix indicating aliased pixels across slices
x      [N 1]      object to be reconstructed
I = S*W*x

I = [I_1; ... I_L]    [L*N 1]       coil images
S = [S_1; ...; S_L]   [L*N MB*N]    coil sensitivity maps   

S_c    [N N]


L  = number of receive coils
MB = multiband factor (number of slices)
N  = number of pixels (n*n, where n = image size)

I      [L*N 1]
S      [L*N MB*N]
```


