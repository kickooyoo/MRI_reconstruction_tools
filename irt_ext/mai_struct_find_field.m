  function [out_name, out_val] = mai_struct_find_field(ss, varargin)
%|function [out_name, out_val] = mai_struct_find_field(ss, varargin)
%|
%| Look recursively through the fields of a structure for a given field.
%| Objects like Fatrix and fatrix2 are converted to structs too for this.
%|
%| Returns ss.(field) or ss.?.(field) or ss.?.?.(field) etc.
%|  
%| Edited by Mai 01/13/2015 for varargin options:
%| - can search by value
%| - can search by value within tolerance 
%| - exhaustive option returns all fields that match query, not just first
%| varargin:
%|		fieldname
%|		name_excludes		cell of strings to exclude in match
%|					e.g. 'image', when searching for 'age'
%|		fieldval
%|		all_vector		default: 1
%|			requires element-by-element match w/i tolerance for each element
%|			if off, returns field if any vector element matches
%|			for use with fieldval
%|		fieldtol		default: 0
%|			for use with fieldval
%|		fieldtol		default: 0
%| 		exhaustive 		default: 0

arg.fieldname = [];
arg.name_excludes = {};
arg.fieldval = [];
arg.fieldtol = 0;
arg.all_vector = 1;
arg.exhaustive = 0;
arg = vararg_pair(arg, varargin);

if nargin == 1 && streq(ss, 'test'), ir_struct_find_field_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

[out_name, out_val, list_name, list_val, ok_name, ok_val] = ir_struct_find_field_do(ss, arg, {}, {});
if ~ok_name && ~ok_val
	fail_msg = 'could not find field';
	if ~isempty(arg.fieldval)
		fail_msg = sprintf('%s with value %d', fail_msg, arg.fieldval);
	end
	if ~isempty(arg.fieldtol)
		fail_msg = sprintf('%s with tolerance %d', fail_msg, arg.fieldtol);
	end
	if ~isempty(arg.fieldname)
		fail_msg = sprintf('%s with name %s', fail_msg, arg.fieldname);
	end
	fail([fail_msg '"%s" in "%s"'], fail_msg, inputname(1));
end
if arg.exhaustive
	out_name = list_name;
	out_val = list_val;
end


function [out_name, out_val, list_name, list_val, ok_name, ok_val] = ir_struct_find_field_do(ss, arg, list_name, list_val)

out_name = [];
out_val = [];
ok_name = isempty(arg.fieldname);
ok_val = isempty(arg.fieldval);
fnames = fieldnames(ss);
for ii=1:length(fnames)
	fname = fnames{ii};

	s1 = ss.(fname);
	if isa(s1, 'Fatrix') || isa(s1, 'fatrix2')
		s1 = struct(s1); % trick
	end
	
	% if match numeric fieldval
	if ~ok_val && ~isempty(s1) && isa(s1, 'numeric') && isa(arg.fieldval, 'numeric')
		out_name = fname;
		out_val = s1;	
		within_tol = abs(double(arg.fieldval(:)) - double(s1(:)))/abs(max(col(arg.fieldval))) <= arg.fieldtol;
		if arg.all_vector
			ok_val = all(within_tol);
		else
			ok_val = any(within_tol);
		end
	end

	% if match str fieldval
	if ~ok_val && ~isempty(s1) && ischar(s1) && ischar(arg.fieldval)
		out_name = fname;
		out_val = s1;	
		ok_val = ~isempty(strfind(lower(s1), lower(arg.fieldval)));
	end
	
	% if match str fieldname
	string_match = ~isempty(strfind(lower(fname), lower(arg.fieldname)));
	if ~isempty(arg.name_excludes)
		for jj = 1:length(arg.name_excludes)
			exclude_match(jj) = ~isempty(strfind(lower(fname), lower(arg.name_excludes{jj})));
		end
		exclude_match = any(exclude_match);
		string_match = string_match && ~exclude_match;
	end
	if ~ok_name && string_match
		out_name = fname;
		out_val = s1;
		ok_name = true;
	end

	if ~arg.exhaustive && ok_val && ok_name
		return;
	elseif ok_val && ok_name
		list_name = [list_name; out_name];
		list_val = [list_val; {out_val}];
		out_name = [];
		out_val = [];
		ok_name = isempty(arg.fieldname);
		ok_val = isempty(arg.fieldval);
	end

	if isstruct(s1) || isa(s1, 'twix_map_obj')
		[out_name, out_val, new_list_name, list_val, ok_name, ok_val] = ir_struct_find_field_do(s1, arg, list_name, list_val);
		for ii = length(list_name)+1:length(new_list_name)
			new_list_name{ii} = sprintf('%s.%s', fname, new_list_name{ii});
		end
		list_name = new_list_name;
		if ok_name && ok_val
			return
		end
	end
end


% ir_struct_find_field_test()
function ir_struct_find_field_test
ss.a = 1;
ss.b.c = 2;
ss.b.d = 3;
jf_equal(1, ir_struct_find_field(ss, 'a'))
jf_equal(3, ir_struct_find_field(ss, 'd'))
try
	ir_struct_find_field(ss, 'bad')
	fail 'should not get here'
catch
end
