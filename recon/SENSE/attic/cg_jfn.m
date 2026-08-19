function [x,res] = cg_jfn(A,b,x_i,nitmax,tol)
% function [x,res] = cg_jfn(A,b,x_i,nitmax,tol)
%
% from CG without the agonizing pain, p 32
%
  
r_i = b - A*x_i;
d_i = r_i;

it = 0;
while norm(r_i) > tol & it < nitmax
	alpha_i = r_i'*r_i/(d_i'*A*d_i);
	x_ii = x_i + alpha_i*d_i;
	r_ii = r_i - alpha_i*A*d_i;
	beta_ii = r_ii'*r_ii/(r_i'*r_i);
	d_ii = r_ii + beta_ii*d_i;

	r_i = r_ii;
	d_i = d_ii;
	x_i = x_ii;
	it = it + 1;

	res(it) = norm(r_i);
end

fprintf(1, 'norm(r) = %.e\n', norm(r_i));
x = x_i;

return
