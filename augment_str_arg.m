function str = augment_str_arg(str, arg, varargin)
% conditionally augment string if certain fields are true
if length(str) >= 4 && strcmp(str(end-3:end), '.mat')
	iarg.ndx_from_end = 4;
else
	iarg.ndx_from_end = 0;
end
iarg.fields_to_add = {'trunc'; 'force_new_sense'; 'whiten'; 'CC'; ...
	'selected_coils'; 'mask_on'; 'tweak_dcf'}; 
iarg = vararg_pair(iarg, varargin);

% vname = @(x) inputname(1);

for ii = 1:length(iarg.fields_to_add)
	if isfield(arg, iarg.fields_to_add{ii});
		curr_field = getfield(arg, iarg.fields_to_add{ii});
		if islogical(curr_field) && curr_field
			str = insert_str(str, iarg.fields_to_add{ii}, iarg.ndx_from_end);
		elseif isnumeric(curr_field) && ~isempty(curr_field)
			str = insert_str(str, iarg.fields_to_add{ii}, iarg.ndx_from_end);
		end
	end
end

end

function new_str = insert_str(str, ins_str, ndx)
% insert ins_str into str at ndx
if ndx >= 1
new_str = [str(1:end-ndx) ins_str str(end-ndx + 1:end)];
else
	new_str = [str ins_str];
end
end

