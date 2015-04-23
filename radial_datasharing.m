function [ds_data, frame_members, ds_freqs, Ns, ds_dcf] = radial_datasharing(freqs, ...
	data, Nyq, varargin)
%function [ds_data, frame_members, ds_freqs, Ns, ds_dcf] = radial_datasharing(freqs, ...
%	data, Nyq, varargin)
%
% generalization of k-Space Weighted Image Contrast (KWIC)
% (essentially data sharing for radial trajectories), but do not assume 
% resample same points of k-space (e.g. Golden Angle)
% if use repeated radial angles, can use conventional datasharing, 
% data_share_fill.m in spline_basis repo
% 
% not yet able to handle multicoil data
%
% inputs:
% freqs: complex-valued, in radians
%		[Nro Nspokespf Nf]
%
% data:	RO values
%		[Nro Nspokes Nf]
%		have to do coil by coil separately!!
% Nyq: maximum azimuthal distance between spokes
%		units: m^-1? radians?
%		leave empty if don't want any datasharing
%
% varargin:
%		Nf:	number of frames
%		Nspokespf: number of spokes per frame
%			need to specify one of the above
%		'format'
%			'logical' default
%			'sparse' TO DO
%			'cells' BUGGY, TO FIX
%		'Fibonnaci' default false, TO DO
%		'Nyquist_spokes', TO DO
%		figs_on
%
% outputs:
% ds_data:
%		{Ns_f Nc}_Nf
%
% frame_members: membership matrices for each frame ndx
%		[Nf Nro_round Nspokes] logical, sparse?
%
% ds_freqs:
%		{Ns_f}_Nf
%		Ns = number of points assigned to frame f after datasharing
%
% Ns: number of samples assigned to each frame, useful for F fatrix
%
% ds_dcf: Voronoi-based density compensation function
%		same size as ds_data
%		
% OR
% frame_members: matrix of cells showing membership
%		[Nro_round Nspokes]
%
% TO DO: Nro_round vs Nro in output!
% To think about: do I use Data at all? should I have a function that lets
% me output datahared data easily and not just member matrices?? UGH
%
% Mai Le, University of Michigan, 01/26/15

if nargin == 1 && streq(freqs, 'test')
	radial_datasharing_test();
	return
elseif nargin == 3 & streq(freqs, 'test')
	radial_datasharing_test(data);
	return
end

% --------------- initializing parameters --------------
% default values for varargin
arg.format = 'logical';
arg.Fibonnaci = false;
arg.Nf = [];
arg.Nspokespf = [];
arg.vary_rings = 0; % as opposed to varying reach
arg.figs_on = 0;
arg.nargout = nargout;
arg = vararg_pair(arg, varargin);
[arg.Nro, arg.Nspokes] = size(freqs);

% check inputs
if nargin < 3, help(mfilename), error(mfilename), end
assert(all(size(freqs) == size(data)), 'freqs and data have mismatched size');
assert(xor(isempty(arg.Nf), isempty(arg.Nspokespf)), ...
	'specify only Nf OR Nspokespf via varargin');
assert((mod(arg.Nf,1) == 0) && (arg.Nf <= arg.Nspokes), ...
	sprintf('invalid Nf: %d', arg.Nf));
assert(~(isempty(arg.Nf) && isempty(arg.Nspokespf)), ...
	'must specify either Nf or Nspokespf');

% calculate respective Nf and Nspokespf
if isempty(arg.Nf)
	arg.Nf = floor(arg.Nspokes/arg.Nspokespf);
end
if isempty(arg.Nspokespf)
	arg.Nspokespf = floor(arg.Nspokes/arg.Nf);
end
if arg.Nf*arg.Nspokespf ~= arg.Nspokes
	display('not yet coded case with leftover spokes!');
	display(sprintf('mismatch: Nspokes (%d) ~= Nspokes/frame (%d) * Nframes (%d)', ...
		arg.Nspokes, arg.Nspokespf, arg.Nf));
	keyboard;
end

% assign data to each original bin
frame_members = trivial_datashare(arg);
if isempty(Nyq)
	display('empty Nyq argument, so no sharing across frames');
	[ds_freqs, ds_data, Ns] = format_outputs(freqs, data, frame_members, arg);
	return; 
else
	init_frame_members = frame_members;
end

% ---------------- datasharing ---------------
% initialize data format
switch arg.format
	case 'logical'
		frame_members = false(arg.Nf, arg.Nro, arg.Nspokes);
	case 'cells'
		frame_members = cell(arg.Nro, arg.Nspokes);
	otherwise
		error(sprintf('unrecognized varargin format %s'), arg.format);
end

% do radial datasharing looping frame by frame
[thetas, data_mags] = Cartesian_to_radial(reshape(freqs, [arg.Nro, arg.Nspokes]));
for frame_ndx = 1:arg.Nf
	
	max_radius = max(col(data_mags));
	
	% get indices of readouts that initially map to this frame
	switch arg.format
		case 'logical'
			frame_theta_ndcs = find(squeeze(init_frame_members(frame_ndx,1,:)) == true);
		case 'cells'
			% simple cells of one val each
			init_frame_members_mat = cell2mat(init_frame_members);
			frame_theta_ndcs = find(init_frame_members_mat(1,:) == frame_ndx);
		otherwise
			error(sprintf('unrecognized varargin format %s'), arg.format);
	end
	
	% do the datasharing
	[ring_thetas, ring_theta_ndcs, annuli] = rdatasharing_1f(thetas, ...
		frame_theta_ndcs, max_radius, Nyq, arg);
	if arg.figs_on
		plot_thetas(ring_thetas, annuli, 'Nyquist', Nyq, 'title', ...
			sprintf('frame %d',frame_ndx));
	end
	
	% format outputs correctly
	switch arg.format
		case 'logical'
			frame_members(frame_ndx,:,:) = format_frame_members(thetas, ...
				data_mags, ring_theta_ndcs, annuli, arg);
		case 'cells'
			tmp = format_frame_members(thetas, data_mags, ring_theta_ndcs, ...
				annuli, arg);
			frame_members = cellfun(@vertcat, frame_members, tmp);
		otherwise
			error(sprintf('unrecognized varargin format %s'), arg.format);
	end	
end

[ds_freqs, ds_data, Ns, ds_dcf] = format_outputs(freqs, data, frame_members, arg);

if (~all(size(ds_dcf) == size(ds_freqs)))
	display('mismatched size for dcf and freqs');
	keyboard;
end

% ----------------- show results -------------
if arg.figs_on
	figure; im(permute(frame_members, [2 3 1]));
end

end

% ??
function [ds_freqs, ds_data, Ns, ds_dcf] = format_outputs(freqs, data, frame_members, arg)
% output columnized freqs and data, remember vary coil last, outside this
% function entirely
	switch arg.format
		case 'logical'
			col_freqs = col(freqs);
			col_data = col(data);
			ds_freqs = [];
			ds_dcf = [];
			ds_data = [];
			Ns = zeros(arg.Nf, 1);
			for frame_ndx = 1:arg.Nf
				curr_members = col(frame_members(frame_ndx,:,:));
				curr_data = col_data(find(curr_members));
				
				if arg.nargout == 5
					curr_freqs = col_freqs(find(curr_members));
					if ~isempty(curr_freqs)
						delta_ro = 1/size(frame_members,2); % normalized freq/Nro
						curr_dcf = calculate_voronoi_dcf(curr_freqs, delta_ro, arg);
					end
					Ns(frame_ndx) = numel(find(curr_members));
					ds_freqs = [ds_freqs; curr_freqs];
					ds_dcf = [ds_dcf; curr_dcf];
				end
				ds_data = [ds_data; curr_data];
			end
		case 'sparse'
			keyboard;
		case 'cells'
			keyboard;
		otherwise
			error(sprintf('unrecognized varargin format %s'), arg.format);
	end
end

% gives voronoi area for a set of frequency points
function dcf = calculate_voronoi_dcf(freqs, delta_ro, arg)
	at_zero = find(freqs == 0);
	uniq_freqs = freqs(setdiff(1:length(freqs),at_zero(2:end)));
	new_zero = find(uniq_freqs == 0);
	if arg.figs_on
		try
			voronoi(real(uniq_freqs), imag(uniq_freqs))
		catch
			keyboard;
		end
	end
	[vor_v, vor_c] = voronoin([real(uniq_freqs) imag(uniq_freqs)]);
	for ii = 1 : size(vor_c ,1)
		ind = vor_c{ii}';
		tess_area(ii,1) = polyarea(vor_v(2:end,1), vor_v(2:end,2));
	end
	if arg.figs_on
		figure;
	end
	for jj = 1:length(vor_c)
		A(jj) = polyarea(vor_v(vor_c{jj},1), vor_v(vor_c{jj},2));
		if isnan(A(jj))
			tmpx = vor_v(vor_c{jj},1);
			tmpy = vor_v(vor_c{jj},2);
			% approximate +- 0.5 boundary
			if length(tmpx) ~= 3
				display('unknown Voronoi shape');
				% hacky
				A(jj) = 0;
			else
			tmpx_trunc = tmpx(~isinf(tmpx));
			tmpy_trunc = tmpy(~isinf(tmpy));
			dists = dist([tmpx_trunc tmpy_trunc]', zeros(2, length(tmpx_trunc)));
			[sorted_dists, dist_ndcs] = sort(dists);
			radius1 = dist([tmpy_trunc(dist_ndcs(1)) tmpx_trunc(dist_ndcs(1))], [0 0]);
			radius2 = dist([tmpy_trunc(dist_ndcs(2)) tmpx_trunc(dist_ndcs(2))], [0 0]);
			new1 = [tmpx_trunc(dist_ndcs(1)) tmpy_trunc(dist_ndcs(1))]*(radius1 + delta_ro)/radius1;
			new2 = [tmpx_trunc(dist_ndcs(2)) tmpy_trunc(dist_ndcs(2))]*(radius2 + delta_ro)/radius2;
			tmpx = cat(1, tmpx_trunc, new1(1), new2(1));
			tmpy = cat(1, tmpy_trunc, new1(2), new2(2));
			if arg.figs_on && (mod(jj,6) == 0)
				scatter(tmpx, tmpy); hold on; 
				scatter(tmpx_trunc, tmpy_trunc); 
				axis([-0.5 0.5 -0.5 0.5])
			end
			if ~all(size(tmpx) == size(tmpy))
				display('size mismatch');
				keyboard;
			end
			A(jj) = polyarea(tmpx, tmpy);
			
			end
		end
	end
	dcf = col(A);
	% add zero value back in
	zero_val = dcf(new_zero);
	dcf(new_zero) = zero_val/numel(at_zero);
	for ii=2:length(at_zero)
		insert = at_zero(ii);
		dcf = [dcf(1:insert-1); dcf(new_zero); dcf(insert:end)];
	end
	if ~all(size(dcf) == size(freqs))
		display('size mismatch with dcf and freqs');
		keyboard;
	end
end

% ??
function frame_members = format_frame_members(thetas, data_mags, ...
	ring_theta_ndcs, radii, arg)
	switch arg.format
		case 'logical'
			ring_members = false(arg.Nro, arg.Nspokes, length(radii));
			for ring_ndx = 1:length(radii)
				correct_spoke = false(arg.Nro, arg.Nspokes);
				correct_spoke(:,ring_theta_ndcs{ring_ndx}) = true;
				correct_annulus = (data_mags <= radii(ring_ndx));
				if ring_ndx > 1
					correct_annulus = correct_annulus & ...
						(data_mags > radii(ring_ndx - 1));
				end
				ring_members(:,:,ring_ndx) = correct_spoke & correct_annulus;
			end
			frame_members = any(ring_members,3);
		case 'sparse'
			keyboard;
		case 'cells'
			ring_members = false(arg.Nro, arg.Nspokespf, length(radii));
			for ring_ndx = 1:length(radii)
				correct_spoke = false(arg.Nro, arg.Nspokes);
				correct_spoke(:,ring_theta_ndcs{ring_ndx}) = true;
				correct_annulus = (data_mags <= radii(ring_ndx));
				if ring_ndx > 1
					correct_annulus = correct_annulus & ...
						(data_mags > radii(ring_ndx - 1));
				end
				ring_members(:,:,ring_ndx) = correct_spoke & correct_annulus;
			end
			ndcs = any(ring_members,3);
			frame_members = cell(arg.Nro, arg.Nspokespf);
			keyboard;
			frame_members(ndcs) =  6;
			
		otherwise
			error(sprintf('unrecognized varargin format %s'), arg.format);
	end
end

% radial datasharing over 1 frame
% combine with previous method?
function [ring_thetas, ring_theta_ndcs, radii] = rdatasharing_1f(thetas, ...
	frame_theta_ndcs, max_radius, Nyq, arg)
	if arg.vary_rings
		[ring_thetas, ring_theta_ndcs, radii] = rdatasharing_1f_set_Nyq( ...
			thetas, frame_theta_ndcs, max_radius, Nyq, arg);
	else
		[ring_thetas, ring_theta_ndcs, radii] = rdatasharing_1f_set_rings(...
			thetas, frame_theta_ndcs, max_radius, Nyq, arg);
	end
end

% datasharing over 1 frame with varying ring size
function [ring_thetas, ring_theta_ndcs, radii] = rdatasharing_1f_set_Nyq(...
	thetas, frame_theta_ndcs, max_radius, Nyq, arg)
	meet_Nyquist = false;
	ring_ndx = 1;
	% add 1 spoke before and after, create new ring
	radii(1) = 0;
	while(radii(end) < max_radius)
		reach = 1;
		while(true)
			if ring_ndx == 1
			ring_theta_ndcs{ring_ndx} = augment_ndx(frame_theta_ndcs, ...
				reach, reach, arg);
			else
			ring_theta_ndcs{ring_ndx} = augment_ndx(ring_theta_ndcs{ring_ndx - 1}, ...
				reach, reach, arg);
			end
			ring_thetas{ring_ndx} = thetas(ring_theta_ndcs{ring_ndx});
			radii(ring_ndx) = Nyquist_radius(ring_thetas{ring_ndx}, Nyq);
			%meet_Nyquist = 
			if (ring_ndx == 1) || (radii(ring_ndx) > radii(ring_ndx - 1))
				break;
			end
			reach = reach + 1;
		end

		if ring_ndx > 1000
			display(sprintf('so many rings?! %d! something buggy...', ...
				ring_ndx));
			keyboard;
		end

		ring_ndx = ring_ndx + 1;
	end
end

% datasharing over 1 frame with set ring size
function [ring_thetas, ring_theta_ndcs, radii] = rdatasharing_1f_set_rings(...
	thetas, frame_theta_ndcs, max_radius, Nyq, arg)
	% set rings to be some preset distance, add spokes as necessary
	% for now, just even sized radii
	meet_Nyquist = false;
	ring_ndx = 1;
	rng(0);
	
	curr_thetas = thetas(frame_theta_ndcs);
	min_radius = Nyquist_radius(curr_thetas, Nyq);
	if min_radius > max_radius
		radii = min_radius;
	else
		radii = min_radius:min_radius:max_radius;
		radii = [radii max_radius];
	end
	Nrings = length(radii);
	% 		ring_theta_ndcs = repmat({frame_theta_ndcs}, 1, Nrings);
	% 		ring_thetas = repmat({curr_thetas}, 1, Nrings);
	for ring_ndx = 1:Nrings
		if ring_ndx == 1
			ring_theta_ndcs{1} = frame_theta_ndcs;
			ring_thetas{1} = curr_thetas;
		else
			ring_theta_ndcs{ring_ndx} = ring_theta_ndcs{ring_ndx - 1};
			ring_thetas{ring_ndx} = ring_thetas{ring_ndx - 1};
		end
		curr_Nyq = Nyquist_radius(ring_thetas{ring_ndx}, Nyq);
		meet_Nyquist = (radii(ring_ndx) <= curr_Nyq);
		counter = 0;
		while(~meet_Nyquist)
			[left_augment, lchange] = augment_ndx(...
				ring_theta_ndcs{ring_ndx}, 1, 0, arg);
			[right_augment, rchange] = augment_ndx(...
				ring_theta_ndcs{ring_ndx}, 0, 1, arg);

			if ~lchange && ~rchange
				display(sprintf('reached all frames on ring %d', ring_ndx));
				break;
			end

			radius_reach_left = azim_dist(thetas(left_augment), radii(ring_ndx));
			radius_reach_right = azim_dist(thetas(right_augment), radii(ring_ndx));
			if (radius_reach_left == radius_reach_right)
				rand_dir = (rand > 0.5);
			else
				rand_dir = NaN;
			end
			if (radius_reach_left < radius_reach_right) || (~isnan(rand_dir) && rand_dir)
				ring_theta_ndcs{ring_ndx} = augment_ndx(...
					ring_theta_ndcs{ring_ndx}, 1, 0, arg);
			elseif (radius_reach_left > radius_reach_right) || (~isnan(rand_dir) && ~rand_dir)
				ring_theta_ndcs{ring_ndx} = augment_ndx(...
					ring_theta_ndcs{ring_ndx}, 0, 1, arg);
			else
				keyboard;
			end
			ring_thetas{ring_ndx} = thetas(ring_theta_ndcs{ring_ndx});
			meet_Nyquist = radii(ring_ndx) <= Nyquist_radius(...
				ring_thetas{ring_ndx}, Nyq);

			counter = counter + 1;
			if counter > arg.Nspokes
				keyboard;
			end
		end
	end
end

% convert Cartesian freqs to radial coordinates
function [thetas, radii] = Cartesian_to_radial(freqs)
% do some averaging to deal with roundoff error
% assume each row of freqs corresponds to one spoke
% NOTE: spoke direction agnostic
	thetas = mod(angle(freqs), pi);
	radii = abs(freqs);
	thetas = mode(thetas,1);
end


function [aug_ndcs, changed] = augment_ndx(ndcs, left, right, arg)
% does auto clipping at [1 arg.Nspokes]
	assert(length(ndcs) > 1, 'ndcs only has one value');
	assert(all(mod([left right],1) == zeros(1,2)), ...
		'invalid augment left/right values');

% check that you only have consecutive indeces
	diffs = ndcs(2:end) - ndcs(1:end-1);
	assert(unique(diffs) == 1, 'current ndcs are nonconsecutive');
	assert(min(ndcs) >= 1, 'invalid lower ndx');
	assert(max(ndcs) <= arg.Nspokes, 'invalid upper ndx');
	
	left_clip = min(left, min(ndcs) - 1);
	right_clip = min(right, arg.Nspokes - max(ndcs));
	aug_ndcs = min(ndcs) - left_clip : max(ndcs) + right_clip;
	
	changed = ~(left_clip == 0 && right_clip == 0);
end

% ??
function Nyq_radius = Nyquist_radius(thetas, Nyq)
% unit agnostic, outputs value in same units as Nyq input
	unit_adist = azim_dist(thetas, 1);
	Nyq_radius = Nyq/unit_adist;
end

% gives azimuthal distance
function max_adist = azim_dist(thetas, radius)
% calculate maximum azimuthal distance between a set of radial spokes of
% fixed radius, assuming that intra-spoke sampling distance is smaller

both_thetas = mod([col(thetas); col(thetas-pi)], 2*pi);
[sorted_thetas, sort_ndcs] = sort(both_thetas);

max_adist = 0;
for gap_ndx = 1:length(thetas)
	a = radius*exp(1i*sorted_thetas(gap_ndx));
	if gap_ndx < length(thetas)
		b = radius*exp(1i*sorted_thetas(gap_ndx + 1));
	else
		b = -radius*exp(1i*sorted_thetas(1));
	end
	curr_dist = dist(a, b);
	if curr_dist > max_adist
		max_adist = curr_dist;
	end
end
end

% Euclidean distance
function dist = dist(x,y)
% good for complex vals
	dist = sqrt(sum(abs(x - y).^2));
end 

% plot annular segments of spokes
function plot_thetas(thetas, radii, varargin)
% assume thetas in cells
% TO DO: enforce same color for spokes in each annulus
% if Nyquist varargin not empty, will plot violations
arg.Nyquist = [];
arg.draw_rings = true;%false;
arg.title = [];

arg = vararg_pair(arg, varargin);

figure;
if length(radii) == 1
	lines = kron([-radii radii], col(exp(1i*thetas{1})));
	plot(lines.');
else
	hold on;
	max_spokes = max(cellfun('size', thetas, 2));
	for ring_ndx = 1:length(radii)
		curr_thetas = thetas{ring_ndx};
		if ring_ndx > 1
			line1 = kron([-radii(ring_ndx) -radii(ring_ndx - 1)], ....
				col(exp(1i*curr_thetas)));
			line2 = kron([radii(ring_ndx - 1) radii(ring_ndx)], ....
				col(exp(1i*curr_thetas)));
		else
			line1 = kron([-radii(ring_ndx) 0], col(exp(1i*curr_thetas)));
			line2 = kron([0 radii(ring_ndx)], col(exp(1i*curr_thetas)));
		end
		plot(line1.');
		hold on;
		plot(line2.');
		if ~isempty(arg.Nyquist)
			both_thetas = mod([col(curr_thetas); col(curr_thetas-pi)], 2*pi);
			[sorted_thetas, sort_ndcs] = sort(both_thetas);
			for gap_ndx = 1:length(curr_thetas)
				a = radii(ring_ndx)*exp(1i*sorted_thetas(gap_ndx));
				if gap_ndx < length(curr_thetas)
					b = radii(ring_ndx)*exp(1i*sorted_thetas(gap_ndx + 1));
				else
					b = -radii(ring_ndx)*exp(1i*sorted_thetas(1));
				end
				curr_dist = dist(a, b);
				if curr_dist - arg.Nyquist > 1e-10
					% tolerance to make up for setting Nyq radius exactly
					% equal
					plot([a; -a; b; -b],'o');
				end
			end
		end
		if arg.draw_rings
			ring = radii(ring_ndx)*exp(1i*linspace(0,2*pi,200));
			plot(ring,'k--');
		end
	end
end
axis equal;
if ~isempty(arg.title) &&ischar(arg.title)
	title(arg.title);
end
end

% trivial case of no datasharing, just mutually exclusive assignment
function frame_members = trivial_datashare(arg)
	switch arg.format
		case {'logical','sparse'}
			in_frame = kron(eye(arg.Nf),ones(arg.Nspokespf,1));
			frame_members = repmat(permute(in_frame, [2 3 1]), [1 arg.Nro 1]);
			if strcmpi(arg.format, 'sparse')
				for frame_ndx = 1:arg.Nf
					frame_members_sparse{frame_ndx} = ...
						sparse(squeeze(frame_members(frame_ndx, :, :)));
				end
				frame_members = frame_members_sparse;
			end
		case 'cells'
			in_frame = ceil((1:arg.Nspokes)/arg.Nspokespf);
			frame_members = repmat(in_frame, [arg.Nro 1]);
			frame_members = mat2cell(frame_members, ones(arg.Nro,1), ones(arg.Nspokes,1));
		otherwise
			error(sprintf('unrecognized varargin format %s'), arg.format);
	end
	
end

% demo/test method
function radial_datasharing_test(varargin)
	arg.datapath = '.';
	arg = vararg_pair(arg, varargin);
	synthetic = 0;
	if synthetic
		arg.Nspokes = 48;
		arg.Nf = 6;
		Nyq = 1.5;
		rng(2);
		thetas = 2*pi*rand(1, arg.Nspokes);

		% Calculate angles for Golden-Angle mode
		GA = 111.246117975/180*pi;
		thetas = [pi/2:GA:GA*arg.Nspokes];
		thetas = mod(thetas,2*pi);
		freqs = kron(col([-10:0.2:10]), exp(1i*thetas));

		[ds_data, frame_members, ds_freqs, Ns, ds_dcf] = radial_datasharing(freqs, ...
			rand(size(freqs)), Nyq, 'Nf', arg.Nf, 'figs_on', 1);
	else
	% GRASP patient data, put datapath in 2nd field (data)
		datafile = [arg.datapath '/XDGRASP_patient1_fast.mat'];
		if exist(datafile, 'file')
			load(datafile);
		end
		
		Nyq = 0.05;
		Nf = 10;
		trunc = 100; % number of spokes, must be <= 1000
		params.Nspokes = trunc;
		
		% do one coil at a time, only do full calc for 1st
		coil_ndx = 1;
		[ds_data(:, coil_ndx), frame_members(:,:,:,coil_ndx), ds_freqs(:, coil_ndx), ...
				 ds_Ns(:, coil_ndx), ds_dcf(:, coil_ndx)] = ...
				radial_datasharing(k(:,1:trunc), ...
				data(:,1:trunc,coil_ndx), Nyq, 'Nf', Nf, 'figs_on', 1);
		for coil_ndx = 2:params.Nc
			[ds_data(:, coil_ndx)] = ...
				radial_datasharing(k(:,1:trunc), ...
				data(:,1:trunc,coil_ndx), Nyq, 'Nf', Nf);
			display(sprintf('done with coil %d/%d', coil_ndx, params.Nc))
		end
		figure; im(permute(frame_members, [2 3 1]));
		% to do: find way to assign data to a frame members of another coil

		F = GsplineF_NC(ds_freqs, ds_Ns, params.Nro, params.Nspokes, Nf, params.Nx, params.Ny, params.Nc);%!!!!
		sos = sos_combine(F'*(repmat(ds_dcf, 1, params.Nc).*ds_data),[],[]);
		figure; im(sos); 
		title('with datasharing and Voronoi dcf');
		figure; im(sos_combine(F'*ds_data,[],[])); 
		title('with datasharing and NO dcf');
	end
end
