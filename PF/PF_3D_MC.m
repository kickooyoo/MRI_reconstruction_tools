function [data, clip_spoke] = PF_3D_MC(data, varargin)
%function [data, params] = PF_3D_MC(data, params)
% wrapper for PF_3D for multicoil
% 
% inputs:
% 	data [Nro Nc Nspokes Nslice_PF]
% varargin:
%	params 
%		struct with lots of info, but specifically Nslice
%	Nslice
%	full_dims: [Nro Nspokes Nslice]
%	block_size (int)
%		break up PF_3D computation by block_size spokes
%		default: 100
%
% outputs:
% 	data [Nro Nc Nspokes_new Nslice]
% 	clip_spoke
%		possibly updated Nspokes (even)
%
% PF_3D option 'PF_location' == [0 0 1] to fit GRASP Siemens data
%	which is radial stack of stars
arg.params = [];
arg.Nslice = [];
arg.full_dims = [];
arg.block_size = 100;
arg.parfor = false;
arg = vararg_pair(arg, varargin);

assert(ndims(data) == 4, 'incorrect data format for PF_3D_MC');
[Nro, Nc, Nspokes, Nslice_PF] = size(data);

if ~isempty(arg.params)
	if isfield(arg.params,'Nslice')
		arg.full_dims = [Nro Nspokes arg.params.Nslice];
	else
		error('invalid params struct is missing field Nslice');
	end
elseif ~isempty(arg.Nslice)
	arg.full_dims = [Nro Nspokes arg.Nslice];
elseif (ndims_ns(arg.full_dims) == 1) && all(sort(size(arg.full_dims)) == [1 3])
	assert(all(arg.full_dims(1:2) == [Nro Nspokes]), ...
		'Nro and Nspokes in full_dims do not match size of data');
else
	display('Invalid set of varargin for PF_3D_MC, choose either:');
	display('(a) params, (b) Nro, Nspokes, Nslice, or (c) full_dims');
	keyboard;
end
	
% enforce Nspokes even
clip_spoke = (mod(arg.full_dims(2), 2) == 1);
data_PF = data(:, :, 1:arg.full_dims(2), :);

data = cat(4, data, zeros(Nro, Nc, Nspokes, arg.full_dims(3) - Nslice_PF));
%full_data = zeros(Nro, Nc, Nspokes, arg.full_dims(3));


all_spoke_ndcs = [];
full_dims = arg.full_dims;

if arg.parfor
	pool = gcp('nocreate');
	if numel(pool) == 0
		pool = parpool();
	end
	parfor coil_ndx = 1:Nc
		coil_data = zeros(Nro, Nspokes, arg.params.Nslice);
		for spoke_block = 1:ceil(Nspokes/arg.block_size)
			if spoke_block == ceil(Nspokes/arg.block_size)
				spoke_ndcs = arg.block_size*(spoke_block - 1) + 1:Nspokes;
			else
				spoke_ndcs = (1:arg.block_size) + arg.block_size*(spoke_block - 1);
			end
			if (mod(numel(spoke_ndcs),2) ~= 0), display('should be even for PF'), keyboard, end
			all_spoke_ndcs = [all_spoke_ndcs spoke_ndcs];
			coil_data_PF = squeeze(data_PF(:, coil_ndx, spoke_ndcs, :));
			curr_dim_2 = numel(spoke_ndcs);
			[~, coil_data(:,spoke_ndcs,:)] = PF_3D(coil_data_PF, [full_dims(1) curr_dim_2 full_dims(3)], ...
				'PF_location', [0 0 1], 'window_step3', 3);
		end
		data(:,coil_ndx, :,:) = coil_data; 
	end
else
	for coil_ndx = 1:Nc
		coil_data = zeros(Nro, Nspokes, arg.params.Nslice);
		for spoke_block = 1:ceil(Nspokes/arg.block_size)
			if spoke_block == ceil(Nspokes/arg.block_size)
				spoke_ndcs = arg.block_size*(spoke_block - 1) + 1:Nspokes;
			else
				spoke_ndcs = (1:arg.block_size) + arg.block_size*(spoke_block - 1);
			end
			if (mod(numel(spoke_ndcs),2) ~= 0), display('should be even for PF'), keyboard, end
			all_spoke_ndcs = [all_spoke_ndcs spoke_ndcs];
			coil_data_PF = squeeze(data_PF(:, coil_ndx, spoke_ndcs, :));
			curr_dim_2 = numel(spoke_ndcs);
			[~, coil_data(:,spoke_ndcs,:)] = PF_3D(coil_data_PF, [full_dims(1) curr_dim_2 full_dims(3)], ...
				'PF_location', [0 0 1], 'window_step3', 3);
		end
		data(:,coil_ndx, :,:) = coil_data; 
	end
end

