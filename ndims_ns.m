function num_non_singleton_dims = ndims_ns(in)
%function num_non_singleton_dims = ndims_ns(in)
% ndims except it returns 1 for vectors, 0 for scalars
num_non_singleton_dims = sum(size(in) > 1);