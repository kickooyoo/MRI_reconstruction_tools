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
%		have to do coil by coil separately!
%
% Nyq: maximum azimuthal distance between spokes
%		units: samples^-1, heuristic: choose 1/max(Nx, Ny)
%		leave empty if don't want any datasharing
%
% varargin:
%		Nf:	number of frames
%		Nspokespf: number of spokes per frame
%			need to specify one of the above
%		'Fibonnaci' enforce center annulus has Fibonnaci number spokes, TO DO
%		'Nyquist_spokes', TO DO
%		'vary_rings' switches between two modes of expanding annuli
%		figs_on
%
% outputs:
% ds_data:
%		[Nds] 
%
% frame_members: membership matrices for each frame ndx
%		[Nf Nro Nspokes] logical
%
% ds_freqs:
%		[Nds], Nds = sum(Ns_f) from f = 1:Nf
%		Ns_f = number of points assigned to frame f after datasharing
%
% Ns:	
%		[Nf] number of samples assigned to each frame, useful for F fatrix
%
% ds_dcf: (optional output)
%		[Nds]
%		Voronoi-based density compensation function
%
% Mai Le, University of Michigan, 01/26/15
% cleaned up 06/30/15

if nargin == 1 && streq(freqs, 'test')
	radial_datasharing_test();
	return
elseif nargin == 3 & streq(freqs, 'test')
	radial_datasharing_test(data);
	return
end

% --------------- initializing parameters --------------
% default values for varargin
arg.Nf = [];
arg.Nspokespf = [];
arg.vary_rings = false; % as opposed to varying reach
arg.figs_on = false;
arg.nargout = nargout;
arg = vararg_pair(arg, varargin);
[arg.Nro, arg.Nspokes] = size(freqs);
arg.Nyq = Nyq;

% check inputs
if nargin < 3, help(mfilename), error(mfilename), end
assert(all(size(freqs) == size(data)), 'freqs and data have mismatched size');
assert(xor(isempty(arg.Nf), isempty(arg.Nspokespf)), ...
	'specify only Nf OR Nspokespf via varargin');
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
% 	display('not yet coded case with leftover spokes!');
% 	display(sprintf('mismatch: Nspokes (%d) ~= Nspokes/frame (%d) * Nframes (%d)', ...
% 		arg.Nspokes, arg.Nspokespf, arg.Nf));
% 	keyboard;
end
assert((mod(arg.Nf,1) == 0) && (arg.Nf <= arg.Nspokes), ...
	sprintf('invalid Nf: %d', arg.Nf));

% assign data to each original bin
frame_members = trivial_datashare(arg);
if isempty(Nyq)
	display('empty Nyq argument, so no sharing across frames');
	[ds_freqs, ds_data, Ns] = format_outputs(freqs, data, frame_members, arg);
	return; 
else
	init_frame_members = frame_members;
end

% --------------------- datasharing ---------------
frame_members = false(arg.Nf, arg.Nro, arg.Nspokes);

% do radial datasharing looping frame by frame
[thetas, data_mags] = Cartesian_to_radial(reshape(freqs, [arg.Nro, arg.Nspokes]));

tic
for frame_ndx = 1:arg.Nf
	
	% determine sample furthest from origin for stopping condition
	arg.max_radius = max(col(data_mags));
	
	% get indices of readouts that initially map to this frame
	frame_theta_ndcs = find(squeeze(init_frame_members(frame_ndx,1,:)) == true);
	
	% do the datasharing
	[ring_thetas, ring_theta_ndcs, annuli] = rdatasharing_1f(thetas, ...
		frame_theta_ndcs, Nyq, arg);
	if arg.figs_on
		plot_thetas(ring_thetas, annuli, 'Nyquist', Nyq, 'title', ...
			sprintf('frame %d',frame_ndx));
	end
	
	% format outputs correctly
	frame_members(frame_ndx,:,:) = format_frame_members(thetas, ...
		data_mags, ring_theta_ndcs, annuli, arg);
end
ds_time = toc;
display(sprintf('done with datasharing in %d sec', ds_time));


% ----------------- show results -------------
if arg.figs_on
	figure; im(permute(frame_members, [2 3 1]));
	title('frame membership');
end

% ----------------- calc dcf, format datashared vector data -------------

[ds_freqs, ds_data, Ns, ds_dcf] = format_outputs(freqs, data, frame_members, arg);

if (arg.nargout == 5) && (~all(size(ds_dcf) == size(ds_freqs)))
	display('mismatched size for dcf and freqs');
	keyboard;
end

end

function [ds_freqs, ds_data, Ns, ds_dcf] = format_outputs(freqs, data, frame_members, arg)
% output columnized freqs and data
	col_freqs = col(freqs);
	col_data = col(data);
	ds_freqs = [];
	ds_dcf = [];
	ds_data = [];
	Ns = zeros(arg.Nf, 1);
	for frame_ndx = 1:arg.Nf
		curr_members = col(frame_members(frame_ndx,:,:));
		curr_data = col_data(find(curr_members));
		curr_freqs = col_freqs(find(curr_members));
		Ns(frame_ndx) = numel(find(curr_members));
		ds_freqs = [ds_freqs; curr_freqs];
		if (arg.nargout == 5) && ~isempty(curr_freqs)
			delta_ro = 1/size(frame_members,2); % normalized freq/Nro
			tic
			curr_dcf = calculate_voronoi_dcf(curr_freqs, delta_ro, arg);
			toc_Voronoi = toc;
			ds_dcf = [ds_dcf; curr_dcf];
			display(sprintf('done with Voronoi for frame %d/%d in % sec', frame_ndx, arg.Nf, toc_Voronoi));
		end
		ds_data = [ds_data; curr_data];
		display(sprintf('done with Voronoi dcf for frame %d/%d', frame_ndx, arg.Nf));
	end
end

% notes to self:
% ring_thetas: 
%	cell array, each cell corresponds to an annulus (indexed
%	center out), values in each cell indicate ANGLE of radial spoke included
%	in each ring, each larger indexed cell should be a superset of any
%	smaller indexed cell
% ring_theta_ndcs:
%	cell array, each cell corresponds to an annulus (indexed
%	center out), values in each cell indicate INDEX of radial spoke included
%	in each ring, each larger indexed cell should be a superset of any
%	smaller indexed cell
% radii:
%	scalar value, radius of each annulus, varies only if arg.vary_rings
function frame_members = format_frame_members(thetas, data_mags, ...
	ring_theta_ndcs, radii, arg)
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
end

% radial datasharing over 1 frame
function [ring_thetas, ring_theta_ndcs, radii] = rdatasharing_1f(thetas, ...
	frame_theta_ndcs, Nyq, arg)
	if arg.vary_rings
		[ring_thetas, ring_theta_ndcs, radii] = rdatasharing_1f_set_Nyq( ...
			thetas, frame_theta_ndcs, Nyq, arg);
	else
		[ring_thetas, ring_theta_ndcs, radii] = rdatasharing_1f_set_rings(...
			thetas, frame_theta_ndcs, Nyq, arg);
	end
end

% datasharing over 1 frame with varying ring size, add 1 spoke before and 
% after, create new ring
function [ring_thetas, ring_theta_ndcs, radii] = rdatasharing_1f_set_Nyq(...
	thetas, frame_theta_ndcs, Nyq, arg)
	meet_Nyquist = false;
	
	% initialize for inner annulus
	ring_theta_ndcs{1} = frame_theta_ndcs';
	ring_thetas{1} = thetas(ring_theta_ndcs{1});
	radii(1) = Nyquist_radius(ring_thetas{1}, Nyq);
	ring_ndx = 2;
	changed = true;
	while (radii(end) < arg.max_radius) && changed % add rings until reach the edge
		first_add = true;
		while(true) % add spokes until reach Nyquist limit within ring
			if first_add
				[new_ndcs, changed] = augment_ndx(ring_theta_ndcs{ring_ndx - 1}, ...
					1, 1, arg);
			else
				[new_ndcs, changed] = augment_ndx(ring_theta_ndcs{ring_ndx}, ...
					1, 1, arg);
			end
			if ~changed
				radii(end) = arg.max_radius;
				break;
			else
				ring_theta_ndcs{ring_ndx} = new_ndcs;
			end
			ring_thetas{ring_ndx} = thetas(ring_theta_ndcs{ring_ndx});
			radii(ring_ndx) = Nyquist_radius(ring_thetas{ring_ndx}, Nyq);
			if (ring_ndx == 1) || (radii(ring_ndx) > radii(ring_ndx - 1))
				break;
			end
			first_add = false;
		end
		if ring_ndx > 100
			display(sprintf('so many rings?! %d! something buggy...', ...
				ring_ndx));
			keyboard;
		end		
		if (length(radii) ~= length(ring_theta_ndcs)) || (length(radii) ~= length(ring_thetas))
			display('mismatched sizes of annuli info');
			keyboard;
		end		
		ring_ndx = ring_ndx + 1;
	end

end

% datasharing over 1 frame with set ring size
function [ring_thetas, ring_theta_ndcs, radii] = rdatasharing_1f_set_rings(...
	thetas, frame_theta_ndcs, Nyq, arg)
	% set rings to be some preset distance, add spokes as necessary
	% for now, just even sized radii
	meet_Nyquist = false;
	ring_ndx = 1;
	rng(0);
	
	curr_thetas = thetas(frame_theta_ndcs);
	min_radius = Nyquist_radius(curr_thetas, Nyq);
	if min_radius > arg.max_radius
		radii = min_radius;
	else
		radii = min_radius:min_radius:arg.max_radius;
		radii = [radii arg.max_radius];
	end
	Nrings = length(radii);
	reached_all_frames = false;
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
				display(sprintf('reached all frames on ring %d/%d', ring_ndx, Nrings));
				reached_all_frames = true;
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
			if counter > 1000; %arg.Nspokes ATTENTION WHY
				keyboard;
			end
		end
		if reached_all_frames && ring_ndx < Nrings
			for xtra_ring_ndx = ring_ndx + 1:Nrings
				ring_theta_ndcs{xtra_ring_ndx} = ring_theta_ndcs{xtra_ring_ndx - 1};
				ring_thetas{xtra_ring_ndx} = ring_thetas{xtra_ring_ndx - 1};
			end
			break;
		end
	end
end

% convert Cartesian freqs to radial coordinates
function [thetas, radii] = Cartesian_to_radial(freqs)
% assume each row of freqs corresponds to one spoke
% NOTE: spoke direction agnostic
	thetas = mod(angle(freqs), pi);
	radii = abs(freqs);
	thetas = mode(thetas,1);
end


function [aug_ndcs, changed] = augment_ndx(ndcs, left, right, arg)
% does auto clipping at [1 arg.Nspokes]
% 	assert(length(ndcs) > 1, 'ndcs only has one value');
	assert(all(mod([left right],1) == zeros(1,2)), ...
		'invalid augment left/right values');

% check that you only have consecutive indeces
	diffs = ndcs(2:end) - ndcs(1:end-1);
	assert(isempty(diffs) || (unique(diffs) == 1), 'current ndcs are nonconsecutive');
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
				if curr_dist >= arg.Nyquist 
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
	trivial_assign = round(linspace(1, arg.Nf, arg.Nspokes))'; 
	in_frame = false(arg.Nspokes, arg.Nf);
	for ii = 1:arg.Nspokes
		in_frame(ii, trivial_assign(ii)) = true;
	end
	frame_members = repmat(permute(in_frame, [2 3 1]), [1 arg.Nro 1]);
end

% demo/test methodF
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
