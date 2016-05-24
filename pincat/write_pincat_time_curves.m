function write_pincat_time_curve(vals, filename)

if length(filename) < 4 || strcmp(filename(end-3:end), '.txt')
	filename = [filename '.txt'];
end

fileID = fopen(filename, 'w');
for ii = 1:length(vals)
	fprintf(fileID, '%.4d\n', vals(ii));
end

fclose(fileID);
