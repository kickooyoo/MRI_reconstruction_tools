function save_im(path, name, varargin)
% function save_im(path, name, varargin)
%
% varargin:
%       FontSize
%       png             default: false

arg.FontSize = [];
arg.png = false;
arg = vararg_pair(arg, varargin);

% if ~isempty(arg.FontSize)
	

print([path name], '-depsc');
savefig([path name])
display(sprintf('saved .fig and .eps named %s in %s', name, path))

if arg.png
        print([path name], '-dpng');
        display(sprintf('saved .png named %s in %s', name, path))

end
