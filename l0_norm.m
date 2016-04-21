function out = l0_norm(in)
% simple l0 norm, no tolerance

out = sum(col(in ~= 0));