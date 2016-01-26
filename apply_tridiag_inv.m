 function output = apply_tridiag_inv(sub, diags, sup, rhs, varargin)
%function output = apply_tridiag_inv(sub, diags, sup, rhs)
%|
%| Solve system of equations with tridiagonal matrix:
%| in:
%|	sub	[n-1]	subdiagonal
%|	diags	[n]	diagonal entries
%|	sup	[n-1]	superdiagonal
%|	rhs	[n, m]	right hand side
%|
%| out
%|	output = T \ rhs, [n, m]
%|
%| 2013-09-15 Mai Le, University of Michigan
%| 2013-09-16 test, comments, white space added by JF

if nargin == 1 && streq(sub, 'test'), ir_apply_tridiag_inv_test, return, end
if nargin < 4, help(mfilename), error(mfilename), end

varg.print = false;
varg = vararg_pair(varg, varargin);

n = length(diags);
assert((length(sub) == n-1) & (length(sup) == n-1) & (size(rhs,1) == n), ...
	'vector lengths incompatible');

% make all inputs col
sub = col(sub);
diags = col(diags);
sup = col(sup);

m = size(rhs,2);

new_sup = zeros(n,1);
new_arg = zeros(n,m);

new_sup(1) = sup(1)/diags(1);
new_arg(1,:) = rhs(1,:)/diags(1);
if varg.print
	ftimer = tic;
end
for ii = 2:n-1
	new_sup(ii) = sup(ii) / (diags(ii) - new_sup(ii-1) * sub(ii-1));
	new_arg(ii,:) = (rhs(ii,:) - new_arg(ii-1,:) * sub(ii-1)) ...
		/ (diags(ii) - new_sup(ii-1) * sub(ii-1));
	if varg.print && mod(ii, floor(n/20)) == 0
		display(sprintf('done with forward pass index %d/%d after %d sec', ii, n, toc(ftimer)));
	end
	for jj = 1:m
		new_arg_2(ii,jj) = (rhs(ii,jj) - new_arg(ii-1,jj) * sub(ii-1)) ...
		/ (diags(ii) - new_sup(ii-1) * sub(ii-1));

	end
end
new_arg(n,:) = (rhs(n,:) - new_arg(n-1,:) * sub(n-1)) ...
	/ (diags(n) - new_sup(n-1) * sub(n-1));

output = zeros(n,m);
output(n,:) = new_arg(n,:);
output_2 = zeros(n,m);
output_2(n,:) = new_arg(n,:);
if varg.print
	rtimer = tic;
end
for ii = n-1:-1:1
	output(ii,:) = new_arg(ii,:) - new_sup(ii) * output(ii+1,:);
	if varg.print && mod(ii, floor(n/1000)) == 0
		display(sprintf('done with reverse pass index %d/%d after %d sec', n-ii, n, toc(rtimer)));
	end
	for jj = 1:m
		output_2(ii,jj) = new_arg(ii,jj) - new_sup(ii) * output_2(ii+1, jj);
	end
end

end



% ir_apply_tridiag_inv_test
% compare to \ as a simple test
function ir_apply_tridiag_inv_test

rng(7)
n = 18;
m = 16;
d = 2 + randn(n, 1);
sub = randn(n-1, 1);
sup = randn(n-1, 1);
rhs = randn(n,m);
tri = diag(d) + diag(sub,-1) + diag(sup,1);
pr cond(tri)
xa = apply_tridiag_inv(sub, d, sup, rhs);
xb = tri \ rhs;
equivs(xa, xb)
%jf_equal(xa, xb)

end

