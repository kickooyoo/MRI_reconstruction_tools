function coutput = arg_struct_to_cell(arg)
%function coutput = arg_struct_to_cell(arg)
% helper function to pass in a struct of varargins
fields = fieldnames(arg);
coutput = {};
for ii = 1:length(fields)
	coutput = [coutput; fields{ii}];
	coutput = [coutput; getfield(arg, fields{ii})];
end

