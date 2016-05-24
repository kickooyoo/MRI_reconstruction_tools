function pincat = load_pincat(dir_name, varargin)
%function pincat = load_pincat(dir_name, varargin)

arg.dims = [256 256 100];
arg.prefix = 'pincat_act';
arg = vararg_pair(arg, varargin);

listing = dir(dir_name);

for ii = 1:length(listing)
	if ~listing(ii).isdir
		curr_file = listing(ii).name;
		ndx = sscanf(curr_file, [arg.prefix '_%d.bin']);
		if ~isempty(ndx) && isempty(strfind(curr_file, 'log'))
			fileID = fopen(curr_file, 'r');
			curr_img = fread(fileID, 'float32');	
			if numel(curr_img) ~= prod(arg.dims), keyboard, end
			curr_img = reshape(curr_img, arg.dims);
			pincat(:,:,:,ndx) = curr_img;
			fclose(fileID);
		end
	end
end

%figure; im('mid3', squeeze(pincat(:,:,:,1)))
