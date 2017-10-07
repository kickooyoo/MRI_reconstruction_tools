function [isFibonacci, varargout] = isFibonacci(in, varargin)
% is the incoming number a Fibonacci number?
% useful in Golden Angle radial sampling
% varargin: 
%	'make_Fibonacci', default: false
arg.make_Fibonacci = false;
arg = vararg_pair(arg, varargin);

list = [1 1];
while list(end) <= in
	list = [list list(end)+list(end-1)];
end
isFibonacci = ismember(in, list);
if arg.make_Fibonacci
	dist = abs(list - in);
	match_ndx = find(dist == min(dist));
	varargout{1} = list(match_ndx);
end

end