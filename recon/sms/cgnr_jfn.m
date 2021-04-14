function [x,res] = cgnr_jfn(A,b,x_i,nitmax,tol)
% function [x,res] = cgnr_jfn(A,b,x_i,nitmax)
%
% Solve Ax = b for rectangular A, using CG on normal equations.
% A must support A*x and A'*y operations.

if ~exist('tol', 'var')
	tol = 1e-6;
end

b = A'*b;   % zero-filled image

r_i = b - AtransAx(A,x_i);
d_i = r_i;

it = 0;
while norm(r_i)/norm(b) > tol & it < nitmax 
	alpha_i = r_i'*r_i/(d_i'*AtransAx(A,d_i));
	x_ii = x_i + alpha_i*d_i;
	r_ii = r_i - alpha_i*AtransAx(A,d_i);
	beta_ii = r_ii'*r_ii/(r_i'*r_i);
	d_ii = r_ii + beta_ii*d_i;

	r_i = r_ii;
	d_i = d_ii;
	x_i = x_ii;
	it = it + 1;

	res(it) = norm(r_i);
end

fprintf(1, '%d CG iterations, norm(r) = %.e\n', it, norm(r_i));

x = x_i;

return;

function AAx = AtransAx(A,x)

Ax = A*x;
AAx = A'*Ax;
% note that A'A is never explicitly formed --> A'A is less sparse than A, which is bad?
% see CG without the pain

return;


