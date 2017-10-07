function out = mod_obi(in, N)
%function out = mod_obi(in, N)
% 
% simple utility for mod for Matlab's one-based indexing
% instead of zero, returns, N

out = mod(in, N);

out(find(out == 0)) = N;
