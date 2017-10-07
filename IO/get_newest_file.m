function fname = get_newest_file(regexp)
% function fname = get_newest(regexp)
% returns filename of the newest file fitting the input regular expression
%
% to do: parse off directories, allow file names to have paths

listing = dir('.');
fname = [];
fdate = 0;
for ii = 1:length(listing)
        if ~listing(ii).isdir
		curr_name = listing(ii).name;
		exp_match = regexpi(curr_name, regexp);
		if ~isempty(exp_match)
			curr_date = listing(ii).datenum;
			if curr_date > fdate
				fname = curr_name;
				fdate = curr_date;
			end
		end
        end
end
