function save_im(path, name, varargin)
% function save_im(path, name, varargin)
%
% varargin:
%       FontSize

arg.FontSize = [];
arg = vararg_pair(arg, varargin);

% if ~isempty(arg.FontSize)
	

print([path name], '-depsc');
savefig([path name])
display(sprintf('saved .fig and .eps named %s in %s', name, path))


