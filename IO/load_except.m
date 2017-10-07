function load_expr = load_except(fname, var)
% function load_except(fname, var)
load_expr = sprintf('load(''%s'', ''-regexp'', ''^(?!%s$).'')', fname, var);
