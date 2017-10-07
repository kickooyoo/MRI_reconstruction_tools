function middle = crop_outside(img, new_dims, varargin)
%function middle = crop_outside(img, new_dims, varargin)
%
% inputs:
%	img (numerical) [dim1, ..., dimN]
%	new_dims (int) [1 ndims]
% output:
%	middle (numerical) [new_dim1, ..., new_dimN]
%		the middle portion of img
% varargin:
%	'odd', (bool) [1 ndims] breaks ties for odd cropping	
%


dims = size(img);
arg.odd = zeros(1, length(dims)); % if uneven cropping, 0 defaults down, 1 defaults up
arg = vararg_pair(arg, varargin);

assert(all(mod(new_dims, 1) == 0), 'new_dims must have integer values');
assert(length(dims) == length(new_dims), 'new dims must mwatch dimensions of img');

indices = '';
for d = 1:length(dims)
	extra = dims(d) - new_dims(d);
	if d > 1
		indices = [indices ', '];
	end
	curr_indices = ':';
	if extra > 1
		if mod(extra, 2) == 0
			curr_indices = sprintf('%d:%d', extra/2+1, dims(d) - extra/2);
		else
			curr_indices = sprintf('%d:%d', (extra + arg.odd(d))/2+1, dims(d) - (extra + arg.odd(d))/2);
		end
	elseif extra < 0
		printf('new dim %d larger than existing size %d, no cropping along dim %d \n', new_dims(d), dims(d), d);
	else % extra = 0, no cropping, no action
	end
	indices = [indices curr_indices];
end

try
	eval(sprintf('middle = img(%s);', indices));
catch
	display('error with eval, check indices');
	keyboard
end

