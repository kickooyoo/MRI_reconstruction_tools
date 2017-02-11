function Nspokes2D = visualize_radial_datasharing(thetas, radii, ...
        frame_members, freqs, varargin)
% thetas: {cell}
% radii: {cell}
% frame_members: membership matrices for each frame ndx
%		[Nf Nro Nspokes] logical
%
arg.Nyquist = [];
arg.draw_rings = true;%false;
arg.title = [];
arg.do_plot = [1 1 1 1];
arg = vararg_pair(arg, varargin);

if nargin == 1 && strcmp(thetas, 'test')
        visualize_radial_datasharing_test();
        return;
end
Nf = length(radii);

Nplots = sum(arg.do_plot);
Msp = floor(sqrt(Nplots));
Nsp = ceil(Nplots/Msp);
for ff = 1:Nf
        figure;
        subplot(Msp, Nsp,1);
        
        if arg.do_plot(1)
                plot_datasharing(thetas{ff}, radii{ff}, arg);
                subplot(Msp, Nsp,2);
        end
        
        if arg.do_plot(2)
                plot_Voronoi(freqs, frame_members, ff);
                subplot(Msp, Nsp,3); % do smarter
        end
        
        if arg.do_plot(3)
                im(permute(frame_members(ff, :,:), [2 3 1]));
                
                subplot(Msp, Nsp,4);
        end
        
        if arg.do_plot(4)
                Nspokes2D(:,:,ff) = plot_temporal_resolution(radii{ff}, frame_members, ff);
        else
                Nspokes_2D = [];e
        end
end

% calculate mean temporal resolution

end

function plot_datasharing(thetas, radii, arg)
% to assign same color at each annulus, need master list of thetas
colors = 'rgbcmyk';
uthetas = [];
for ring_ndx = 1:length(thetas)
        uthetas = unique([uthetas; col(thetas{ring_ndx})]);
end

if length(radii) == 1
        lines = kron([-radii radii], col(exp(1i*thetas{1})));
        plot(lines.');
else
        hold on;
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
                hold on;
                color_order = vfind(curr_thetas, uthetas);
                cmplot(line1.', colors, color_order);
                cmplot(line2.', colors, color_order);
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
axis tight;
if ~isempty(arg.title) &&ischar(arg.title)
        title(arg.title);
end

end

% plots amatrix as lines with specific color order
function cmplot(mat, colors, color_order)
[M, N] = size(mat);
for ii = 1:N
        curr_color = colors(mod(color_order(ii) - 1, length(colors)) + 1);
        plot(mat(:,ii), curr_color);
end
end

function indices = vfind(query, database)

for ii = 1:length(query)
        indices(ii) = find(database == query(ii));
end
end


function plot_Voronoi(freqs, frame_members, ff)
% is this monkey business necessary?
col_freqs = col(freqs);
curr_members = col(frame_members(ff,:,:));
curr_freqs = col_freqs(find(curr_members));
first_instance = zeros(1,length(curr_freqs));
for ii = 1:length(curr_freqs)
        all_instances = find(curr_freqs == curr_freqs(ii));
        if length(all_instances) < 1, keyboard; end
        if max(all_instances) > length(curr_freqs), keyboard; end
        first_instance(ii) = all_instances(1);
end
uniq_freqs = curr_freqs(find(first_instance == (1:length(curr_freqs))));
uniq_angles = unique(angle([uniq_freqs; -uniq_freqs]));
max_radius = max(col(abs(curr_freqs)));
delta_ro = 1/size(frame_members,2);
ring = unique((max_radius + delta_ro)*exp(1*i*uniq_angles));
uniq_freqs = [uniq_freqs; ring];
voronoi(real(uniq_freqs), imag(uniq_freqs))
axis equal;
axis tight;

end

function Nspokes2D = plot_temporal_resolution(radii, frame_members, ff)
a = linspace(- max(radii), max(radii), 500);
[xx, yy] = ndgrid(a, a);
Nro = size(frame_members, 2);
ro = [-Nro/2:Nro/2-1]/Nro;
dist2D = sqrt(xx.^2 + yy.^2);
Nspokes = squeeze(sum(frame_members(ff, :, :), 3));
Nspokes2D = interp1(ro, Nspokes, dist2D, 'nearest' , 'extrap');
im(Nspokes2D);

end

% Euclidean distance
function dist = dist(x,y)
dist = sqrt(sum(abs(x - y).^2));
end


function  visualize_radial_datasharing_test()

Nf = 10;
Nspokes = 20; % per frame
Nro = 50;
[k, phi] = get_GA_coords(Nro, Nspokes*Nf, 1);
k = reshape(k, Nro, Nspokes, Nf);
data = randn(Nro, Nspokes, Nf);
Nyq = 0.05;
[ds_data, frame_members, ds_freqs, Ns, ds_dcf] = radial_datasharing(k, ...
        data, Nyq, Nf, 'figs_on', 1);
% figure; jf_slicer_ml(frame_members)
% visualize_radial_datasharing(frame_members);

keyboard

end