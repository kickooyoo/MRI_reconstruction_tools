function cov_mat = calc_coil_cov(noise, varargin);
%function cov_mat = calc_coil_cov(noise, varargin);
arg.opt = 'coeff';% opt for xcov: 'biased', 'unbiased', 'coeff', 'none'
arg = vararg_pair(arg, varargin);

[Nc, Nmeas] = size(noise);

for ii = 1:Nc
	for jj = 1:Nc
		cov_mat(ii, jj) = xcov(noise(ii,:), noise(jj,:), 0, arg.opt);
	end
end
