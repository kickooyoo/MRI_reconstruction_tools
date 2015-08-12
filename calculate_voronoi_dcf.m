function dcf = calculate_voronoi_dcf(freqs, delta_ro, arg)
%function dcf = calculate_voronoi_dcf(freqs, delta_ro, arg)
%
% gives voronoi area for a set of frequency points, allows for duplicate
% points (e.g. at DC)
%
% pad ring makes an evenly spaced ring but can create weird edge polygons
arg.pad_ring = false;
arg.pad_spokes = true;
arg.check_polygon = false;

% pare down all frequencies into set of unique points
% cannot use Matlab's unique because need to reexpand at end
for ii = 1:length(freqs)
	all_instances = find(freqs == freqs(ii));
	if length(all_instances) < 1, keyboard; end
	if max(all_instances) > length(freqs), keyboard; end
	first_instance(ii) = all_instances(1);
end
uniq_freqs = freqs(find(first_instance == (1:length(freqs))));
Nuniq = length(uniq_freqs);

if arg.pad_ring
	max_radius = max(abs(col(uniq_freqs)));
	Npad = 4*sum(abs(abs(freqs) - arg.max_radius) < delta_ro*1e-4);
	ring = (arg.max_radius + delta_ro)*exp(1*i*linspace(0,2*pi,Npad)');
	uniq_freqs = [uniq_freqs; ring];
elseif arg.pad_spokes
	uniq_angles = unique(angle([uniq_freqs; -uniq_freqs]));
	Npad = numel(uniq_angles);
	ring = (arg.max_radius + delta_ro)*exp(1*i*uniq_angles);
	uniq_freqs = [uniq_freqs; ring];
end

% use Matlab's built in Voronoi commands to find area of nearest neighbor
% cells
if arg.figs_on
	try
		figure;
		voronoi(real(uniq_freqs), imag(uniq_freqs))
		axis equal;
		axis tight;
	catch
		display('error when attempting to plot voronoi cells');
		keyboard;
	end
end
[vor_v, vor_c] = voronoin([real(uniq_freqs) imag(uniq_freqs)]);

% for each centroid, calculate area from vertices of polygon
% only go up to Nuniq, remainder points are for padding ring
for jj = 1:Nuniq 
	A(jj) = polyarea(vor_v(vor_c{jj},1), vor_v(vor_c{jj},2));
	tmpx = vor_v(vor_c{jj},1);
	tmpy = vor_v(vor_c{jj},2);
	if arg.pad_ring
		% jagged polygons often go over boundary
		boundary = arg.max_radius + 10*delta_ro + arg.Nyq/10; 
	elseif arg.pad_spokes
		boundary = arg.max_radius + delta_ro + arg.Nyq/10;
	else
		boundary = arg.max_radius;
	end
	out_of_bounds = any(dist([tmpx'; tmpy'], zeros(2, size(tmpx,1))) > boundary);
	if (isnan(A(jj)) || out_of_bounds) && arg.check_polygon
		display('invalid polygon');
		keyboard;
	end
	if (mod(jj,1000) == 0), display(sprintf('done with Voronoi areas for %d/%d', jj, Nuniq)); end
end
uniq_dcf = col(A);

uniq_spread_dcf = embed(uniq_dcf, first_instance == (1:length(freqs)));
for ii = 1:length(freqs)
	if first_instance(ii) > length(uniq_spread_dcf), keyboard; end
	num_occurences = length(find(first_instance == first_instance(ii)));
	dcf(ii) = uniq_spread_dcf(first_instance(ii))/num_occurences;
end
dcf = col(dcf);

if ~all(size(dcf) == size(freqs))
	display('size mismatch with dcf and freqs');
	keyboard;
end
end

% Euclidean distance
function dist = dist(x,y)
% good for complex vals
dist = sqrt(sum(abs(x - y).^2));
end