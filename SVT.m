function X = SVT(M)
% implementation of SVT by CAndes

% assume M fully sampled
X = M;

[N1, N2] = size(X)
Y = zeros(N1, N2);

counter = 0;
max_iter = 20;
stop_condition = counter > max_iter;
step = ones(1,max_iter);

while(~stop_condition)
	counter = counter + 1;
	X = soft_thresh_mat(Y);
	Y = Y + step(counter)*orth_proj(M-X);
	
end
end


function out = soft_thresh_mat(in)
[U, S, V] = svd(in);
thresh = 0.005*S(1,1);
S_soft_thresh = max(S - thresh,0);
out = U*S_soft_thresh*V';
end

function out = orth_proj(in)

end