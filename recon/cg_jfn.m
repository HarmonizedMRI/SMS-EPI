% from CG without the agonizing pain, p 32
%
% $Id: cg_jfn.m,v 1.4 2010-09-23 02:32:47 jfnielse Exp $

function [x] = cg_jfn(A,b,x_i,tol)
  
r_i = b - A*x_i;
d_i = r_i;

it = 0;
nit = 10;
while norm(r_i) > tol & it < nit
	alpha_i = r_i'*r_i/(d_i'*A*d_i);
	x_ii = x_i + alpha_i*d_i;
	r_ii = r_i - alpha_i*A*d_i;
	beta_ii = r_ii'*r_ii/(r_i'*r_i);
	d_ii = r_ii + beta_ii*d_i;

	r_i = r_ii;
	d_i = d_ii;
	x_i = x_ii;
	it = it + 1;
end

fprintf(1, 'norm(r) = %.e\n', norm(r_i));
x = x_i;

return
