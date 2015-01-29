function [frame_members, ds_freqs, ds_data, Ns] = radial_datasharing(freqs, data, Nyq, varargin)
%function [frame_members, ds_freqs, ds_data] = radial_datasharing(freqs, data, Nyq, varargin)
% generalization of k-Space Weighted Image Contrast (KWIC)
% (essentially data sharing for radial trajectories), but do not assume 
% resample same points of k-space (e.g. Golden Angle)
% if use repeated radial angles, can use conventional datasharing, 
% data_share_fill.m in spline_basis repo
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
%
% outputs:
% frame_members: membership matrices for each frame ndx
%		[Nf Nro_round Nspokes] logical, sparse?
%
% ds_freqs:
%		{Ns_f}_Nf
%		Ns = number of points assigned to frame f after datasharing
%
% ds_data:
%		{Ns_f Nc}_Nf
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

if nargin == 1 && streq(freqs, 'test'), radial_datasharing_test, return, end

arg.format = 'logical';
arg.Fibonnaci = false;
arg.Nf = [];
arg.Nspokespf = [];
arg.vary_rings = 0; % as opposed to varying reach
arg.figs_on = 0;
arg = vararg_pair(arg, varargin);
[arg.Nro, arg.Nspokes] = size(freqs);

if nargin < 3, help(mfilename), error(mfilename), end
assert(all(size(freqs) == size(data)), 'freqs and data have mismatched size');
assert(xor(isempty(arg.Nf), isempty(arg.Nspokespf)), ...
	'specify only Nf OR Nspokespf via varargin');
assert((mod(arg.Nf,1) == 0) && (arg.Nf <= arg.Nspokes), ...
	sprintf('invalid Nf: %d', arg.Nf));

% calculate respective Nf and Nspokespf
if isempty(arg.Nf)
	arg.Nf = floor(arg.Nspokes/arg.Nspokespf);
end
if isempty(arg.Nspokespf)
	arg.Nspokespf = floor(arg.Nspokes/arg.Nf);
end
if arg.Nf*arg.Nspokespf ~= arg.Nspokes
	display('not yet coded case with leftover spokes!');
	keyboard;
end

frame_members = trivial_datashare(arg);
if isempty(Nyq)
	[ds_freqs, ds_data, Ns] = format_outputs(freqs, data, frame_members, arg);
	return; 
else
	init_frame_members = frame_members;
end

switch arg.format
	case 'logical'
		frame_members = false(arg.Nf, arg.Nro, arg.Nspokes);
	case 'sparse'
		keyboard;
	case 'cells'
		frame_members = cell(arg.Nro, arg.Nspokes);
	otherwise
		error(sprintf('unrecognized varargin format %s'), arg.format);
end

% do radial datasharing frame by frame
[thetas, data_mags] = Cartesian_to_radial(reshape(freqs, [arg.Nro, arg.Nspokes]));
for frame_ndx = 1:arg.Nf
	max_radius = max(col(data_mags));
	switch arg.format
		case 'logical'
			frame_theta_ndcs = find(squeeze(init_frame_members(frame_ndx,1,:)) == true);
		case 'sparse'
			keyboard;
		case 'cells'
			% simple cells of one val each
			init_frame_members_mat = cell2mat(init_frame_members);
			frame_theta_ndcs = find(init_frame_members_mat(1,:) == frame_ndx);
		otherwise
			error(sprintf('unrecognized varargin format %s'), arg.format);
	end
	
	[ring_thetas, ring_theta_ndcs, annuli] = rdatasharing_1f(thetas, ...
		frame_theta_ndcs, max_radius, Nyq, arg);
	if arg.figs_on
		plot_thetas(ring_thetas, annuli, 'Nyquist', Nyq, 'title', ...
			sprintf('frame %d',frame_ndx));
	end
	
	switch arg.format
		case 'logical'
			frame_members(frame_ndx,:,:) = format_frame_members(thetas, ...
				data_mags, ring_theta_ndcs, annuli, arg);
		case 'sparse'
			keyboard;
		case 'cells'
			tmp = format_frame_members(thetas, data_mags, ring_theta_ndcs, ...
				annuli, arg);
			frame_members = cellfun(@vertcat, frame_members, tmp);
		otherwise
			error(sprintf('unrecognized varargin format %s'), arg.format);
	end
	[ds_freqs, ds_data, Ns] = format_outputs(freqs, data, frame_members, arg);
end

if arg.figs_on
	figure; im(permute(frame_members, [2 3 1]));
end

end

function [ds_freqs, ds_data, Ns] = format_outputs(freqs, data, frame_members, arg)
% output columnized freqs and data, remember vary coil last, outside this
% function entirely
	switch arg.format
		case 'logical'
			col_freqs = col(freqs);
			col_data = col(data);
			ds_freqs = [];
			ds_data = [];
			Ns = zeros(arg.Nf, 1);
			for frame_ndx = 1:arg.Nf
				curr_members = col(frame_members(frame_ndx,:,:));
				curr_freqs = col_freqs(find(curr_members));
				curr_data = col_data(find(curr_members));
				Ns(frame_ndx) = numel(find(curr_members));
				ds_freqs = [ds_freqs; curr_freqs];
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

function [ring_thetas, ring_theta_ndcs, radii] = rdatasharing_1f(thetas, ...
	frame_theta_ndcs, max_radius, Nyq, arg)
	meet_Nyquist = false;
	ring_ndx = 1;
	if arg.vary_rings
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
		
	else
		% set rings to be some preset distance, add spokes as necessary
		% for now, just even sized radii
		curr_thetas = thetas(frame_theta_ndcs);
		min_radius = Nyquist_radius(curr_thetas, Nyq);
		if min_radius > max_radius
			radii = min_radius;
		else
			radii = min_radius:min_radius:max_radius;
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
end

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

% check that you only have consectuive indeces
	diffs = ndcs(2:end) - ndcs(1:end-1);
	assert(unique(diffs) == 1, 'current ndcs are nonconsecutive');
	assert(min(ndcs) >= 1, 'invalid lower ndx');
	assert(max(ndcs) <= arg.Nspokes, 'invalid upper ndx');
	
	left_clip = min(left, min(ndcs) - 1);
	right_clip = min(right, arg.Nspokes - max(ndcs));
	aug_ndcs = min(ndcs) - left_clip : max(ndcs) + right_clip;
	
	changed = ~(left_clip == 0 && right_clip == 0);
end

function Nyq_radius = Nyquist_radius(thetas, Nyq)
% unit agnostic, outputs value in same units as Nyq input
	unit_adist = azim_dist(thetas, 1);
	Nyq_radius = Nyq/unit_adist;
end

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

function dist = dist(x,y)
% good for complex vals
	dist = sqrt(sum(abs(x - y).^2));
end 

function plot_thetas(thetas, radii, varargin)
% plot annular segments of spokes
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

function frame_members = trivial_datashare(arg)
% trivial case of no datasharing
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

function radial_datasharing_test()
	%close all;
	arg.Nspokes = 36;
	arg.Nf = 6;
	Nyq = 3;
	rng(2);
	thetas = 2*pi*rand(1, arg.Nspokes);
	
	% Calculate angles for Golden-Angle mode
	GA = 111.246117975/180*pi;
	thetas = [pi/2:GA:GA*arg.Nspokes];
	thetas = mod(thetas,2*pi);
	%figure; plot_thetas({thetas}, 1);
	
	freqs = kron(col([-10:0.2:10]), exp(1i*thetas));
	%figure; scatter(real(freqs(:)), imag(freqs(:)));
	%frame_members = radial_datasharing(freqs, rand(size(freqs)), Nyq, 'Nf', ...
%		arg.Nf, 'format','cells');
	
	[frame_members, ds_freqs, ds_data, Ns] = radial_datasharing(freqs, ...
		rand(size(freqs)), Nyq, 'Nf', arg.Nf);
	keyboard;
	%frame_members = radial_datasharing(freqs, rand(size(freqs)), Nyq, 'Nf', ...
	%	arg.Nf, 'vary_rings', 1);
	%[thetas, radii] = Cartesian_to_radial(freqs);
end
