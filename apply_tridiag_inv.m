function output = apply_tridiag_inv(sub, diags, sup, arg, varargin)
% function output = apply_tridiag_inv(sub, diags, sup, arg, varargin)
% sub: subdiagonal, length n-1
% diags: diagonal entries, length n
% sup: superdiagonal, length n-1
varg.print = false;
varg = vararg_pair(varg, varargin);

n = length(diags);
assert((length(sub) == n-1) & (length(sup) == n-1),'vectors incompatible lengths');

new_sup = zeros(n,1);
new_arg = zeros(n,1);

new_sup(1) = sup(1)/diags(1);
new_arg(1) = arg(1)/diags(1);

if varg.print
	ftimer = tic;
end

for ii = 2:n-1
	new_sup(ii) = sup(ii)/(diags(ii)-new_sup(ii-1)*sub(ii-1));
	new_arg(ii) = (arg(ii)-new_arg(ii-1)*sub(ii-1))/(diags(ii)-new_sup(ii-1)*sub(ii-1));
	if varg.print && mod(ii, floor(n/20)) == 0
		display(sprintf('done with forward pass index %d/%d after %d sec', ii, n, toc(ftimer)));
	end
end
new_arg(n) = (arg(n)-new_arg(n-1)*sub(n-1))/(diags(n)-new_sup(n-1)*sub(n-1));

output = zeros(n,1);
output(n) = new_arg(n);

if varg.print
	rtimer = tic;
end

for ii = n-1:-1:1
	output(ii) = new_arg(ii)-new_sup(ii)*output(ii+1);
	if varg.print && mod(ii, floor(n/1000)) == 0
		display(sprintf('done with reverse pass index %d/%d after %d sec', n-ii, n, toc(rtimer)));
	end
% 	display(sprintf('%d', ii));
end

