function [varg_trim, leftovers] = trim_vararg(varg, arg)

% to use with vararg_pair
% cuts extraneous fields to avoid error with vararg_pair

assert(mod(length(varg),2) == 0, 'must have even number entries corresponding to name-value pairs');

varg_trim = {};
leftovers = {};
for ii = 1:length(varg)/2
	field = varg{1 + (ii - 1)*2};
	value = varg{ii*2};
	assert(ischar(field), sprintf('field %d, %d th entry of varg not string', field, ii));
	if isfield(arg, field)
		varg_trim = [varg_trim; field; value];
	else
		leftovers = [leftovers; field; value];
	end
end



