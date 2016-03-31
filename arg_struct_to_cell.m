function coutput = arg_struct_to_cell(arg)
%function coutput = arg_struct_to_cell(arg)
% helper function to pass in a struct of varargins
fields = fieldnames(arg);
coutput = {};
for ii = 1:length(fields)
	coutput{(ii - 1)*2 + 1} = fields{ii};
	coutput{2*ii} = getfield(arg, fields{ii});
end

