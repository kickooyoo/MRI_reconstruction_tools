function [output, filenames] = exist_regexp(name_regex, varargin)
% 
% similar to Matlab's exist but with ability to detect for specific regex
% helpful for prefixes
% caution: may return too many valid files because matlab regexp behavior
% only tested on prefixes, e.g., 'testing*"

if nargin == 2
        type = varargin{1};
else 
        type = 'any';
end

dir_ndx = strfind(name_regex, '/');
if ~isempty(dir_ndx)
	dir_ndx = dir_ndx(end);
        curr_dir = name_regex(1:dir_ndx-1);
        fname = name_regex(dir_ndx+1:end);
else
        curr_dir = '.';
        fname = name_regex;
end

listing = dir(curr_dir);

output = [];
filenames = {};
counter = 1;
for ii = 1:length(listing)
        
        % to do: file type checking
        
        name_match = regexp(listing(ii).name, fname);
        
        if ~isempty(name_match)
                full_name = [curr_dir '/' listing(ii).name];
                if strcmp(type, 'any')
                        output(counter) = exist(full_name);
                else
                        output(counter) = exist(full_name, type);
                end
                filenames{counter} = full_name;
                counter = counter + 1;
        end
       
        
end
