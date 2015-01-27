function x = ir_wls_init_scale_mai(A, y, varargin)
% function x = ir_wls_init_scale_mai(A, y, varargin)
% varargin x, default is A'*y
arg.x = A' * y;
arg = vararg_pair(arg, varargin);
tmp = A * arg.x;
scale = sum(col(conj(tmp) .* y), 'double') / sum(col(abs(tmp).^2), 'double');
x = scale * arg.x;
