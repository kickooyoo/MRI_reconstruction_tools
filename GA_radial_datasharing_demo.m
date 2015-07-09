% script to demo golden angle k-space data sharing
% 
% calls create_GA_radial_samples.m and radial_datasharing.m 
% but not calculate_voronoi_dcf.m

%% Toy values
Nf = 3;%5;%3;
Nro = 20;
Nx = 15;
Fibo = 8;
Nspokes = Nf*Fibo;
grad_shift = -1.75;

%% somewhat more realistic set of values
Nf = 5;
Nro = 100;
Nx = 64; % assume square image for simplicity
Fibo = 34;
Nspokes = Nf*Fibo; % multiple of Fibonacci number
grad_shift = -0.75; 

%%

k = create_GA_radial_samples(Nspokes, Nro, 'grad_shift', grad_shift, 'figs_on', true);
% figure; plot(k)

% radial_datasharing requires input of data, so using nonsensedata
data = rand(size(k));
Nyq = 1/Nx; 
[ds_data, frame_members, ds_freqs, Ns, ds_dcf] = radial_datasharing(k, ...
	data, Nyq, 'Nspokespf', Fibo, 'figs_on', true, 'vary_rings', true);
% alternatively, [ds_data, frame_members, ds_freqs, Ns] = radial_datasharing(k, ...
%	data, Nyq, 'Nf', Nf, 'figs_on', true);
% note: blue circles indicate spokes where azimuthal distance is within
% 1e-10 of Nyquist parameter
% plotted spokes are simplified and do not show gradient shift effect


%% to display neatly which datashared freqs are assigned to which frame

% frame_members = true(Nf, Nro, Nspokes);
colors = 'rgbkc';
figure; hold on;
for ii = 1:Nf
	% sum(col(frame_members)) == numel(ds_freqs), use indeces of
	% frame_members to extract relevant values for each frame
	frame_mask = false(Nf, Nro, Nspokes);
	frame_mask(ii, :, :) = true(Nro, Nspokes);
	frame_mask = frame_mask & frame_members;
	frame_mask_frame = frame_mask(find(col(frame_members)));
	ds_freqs_frame{ii} = ds_freqs(find(frame_mask_frame));
	hold on; 
	plot(ds_freqs_frame{ii}, sprintf('%s.',colors(ii))); 
	title(sprintf('frame: %d', ii));
	axis equal;
end
