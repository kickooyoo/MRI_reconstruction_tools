function cov_mat = calc_coil_cov(noise)
%function cov_mat = calc_coil_cov(noise)

[Nc, Nmeas] = size(noise);

for ii = 1:Nc
	for jj = 1:Nc
		cov_mat(ii, jj) = xcov(noise(ii,:), noise(jj,:), 0, 'coeff');
	end
end
