function [Nc_sel, varargout] = select_coils(first_coil_dim, selected_coils, varargin)
% e.g. [params.Nc, data, noise_data, sense_maps] = select_coils(data,
% noise_data, sense_maps);
%

assert(nargin >= 3, 'select_coils requires at least one input to truncate');
assert(ndims_ns(size(varargin{1}) >= first_coil_dim, 'first_coil_ndx out of bounds');
Nc = size(varargin{1}, 1);

if any(selected_coils > Nc) || any(selected_coils < 1)
	display('invalied selected coils');
	keyboard;
end

for ii = 1:nargin-2
	curr = varargin{ii};
	coil_dim = find(size(curr) == Nc);
	if ~unique(coil_dim)
		display('non unique inferred coil dim');
		keyboard;
	end
	
	if ndims_ns(curr) > 3
		display('not yet coded for 4D inputs');
		keyboard;
	end
	
	% smarter implemenation with reshape and index counting?
	switch coil_dim
		case 1
			varargout{ii} = curr(selected_coils, :, :);
		case 2
			varargout{ii} = curr(:, selected_coils, :);
		case 3
			varargout{ii} = curr(:, :, selected_coils);
		otherwise
	end
% 	varargout{ii} = curr;
end

Nc_sel = numel(selected_coils);