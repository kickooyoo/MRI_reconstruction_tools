function time_series_clicker(imgs, varargin)
% expect imgs to be dynamic 2D [Nx, Ny, Nt, NI]
arg.texts = {};
arg.texts_x = 10;
arg.texts_y = 10;
arg.colors = {'k--', 'b--', 'r--', 'g--'};
arg.avg_ball = []; % TODO
arg.clim = [0 max(abs(col(imgs)))];
arg = vararg_pair(arg, varargin);

if ischar(imgs) && strcmpi(imgs, 'test')
	time_series_clicker_test();
	return;
end

[Nx, Ny, Nt, NI] = size(imgs);
Nrows = ceil(sqrt(NI + 1));
Ncols = ceil(NI/Nrows);

curr_x = 1;
curr_y = 1;
curr_t = 1;

% start with just first frame
figure;
for ii = 1:NI
	sp(ii) = subplot(Nrows, Ncols, ii);
	im(imgs(:,:, curr_t,ii), arg.clim);
	if ~isempty(arg.texts)
		text(mod(arg.texts_x, Nx), mod(arg.texts_y, Ny), ...
			arg.texts{ii}, 'color', [1 1 1]);
	end
	hold on;
	scatter(curr_x, curr_y, 'r*');
	hold off;
end
sp(NI + 1) = subplot(Nrows, Ncols, NI+1);
hold on;
for ii = 1:NI
	plot(1:Nt, squeeze(imgs(curr_x, curr_y, :, ii)), ...
		arg.colors{mod_obi(ii, length(arg.colors))});
end
% separate loop for legend indexing
for ii = 1:NI
	scatter(curr_t, squeeze(imgs(curr_x, curr_y, curr_t, ii)), 'mo');
end
plot(curr_t*ones(1,100), linspace(min(ylim), max(ylim), 100), 'm--');
axis([1 Nt min(0, min(ylim)) max(ylim)]);
hold off;
title(sprintf('time series at voxel (%d, %d)', curr_x, curr_y));
xlabel('frame index');
ylabel('voxel value');
if ~isempty(arg.texts)
	hlegend = legend(arg.texts);
	set(hlegend, 'Location', 'southoutside', 'FontSize', 8);
end

display('to quit selecting points, press "q"');
while (true)
	
	% to quit, press 'q'
	kkey = get(gcf,'CurrentCharacter');
	if strcmp(kkey, 'q')
		break;
	end
	
	% hack for making crosshairs of ginput visible on Mac
	subplot(Nrows, Ncols, 1);
	hack = legend('');
	set(hack, 'position', [-5 -5 5 5]);

	% get click and active subplot
	[xcoord, ycoord] = ginput(1); 
	clicked_plot = gca;
	active_sp = find(sp == clicked_plot);
	
	if active_sp <= NI % update time series
		curr_x = mod_obi(round(xcoord), Nx);
		curr_y = mod_obi(round(ycoord), Ny);
	else % update frame
		curr_t = min(max(round(xcoord), 1), Nt);
	end
	
	% redraw each pic
	for ii = 1:NI
		sp(ii) = subplot(Nrows, Ncols, ii);
		im(imgs(:,:,curr_t,ii), arg.clim);
		if ~isempty(arg.texts)
			text(mod(arg.texts_x, Nx), mod(arg.texts_y, Ny), ...
				arg.texts{ii}, 'color', [1 1 1]);
		end
		hold on;
		scatter(curr_x, curr_y, 'r*');
		hold off;
	end
	
	% redraw time series
	cla(sp(NI + 1), 'reset');
	sp(NI + 1) = subplot(Nrows, Ncols, NI + 1);
	hold on;
	for ii = 1:NI
		plot(1:Nt, squeeze(imgs(curr_x, curr_y, :, ii)), ...
			arg.colors{mod_obi(ii, length(arg.colors))});
	end
	% separate loop for legend indexing
	for ii = 1:NI
		scatter(curr_t, squeeze(imgs(curr_x, curr_y, curr_t, ii)), 'mo');
	end
	plot(curr_t*ones(1,100), linspace(min(ylim), max(ylim), 100), 'm--');
	axis([1 Nt min(0, min(ylim)) max(ylim)]);
	title(sprintf('time series at voxel (%d, %d)', curr_x, curr_y));
	xlabel('frame index');
	ylabel('voxel value');
	if ~isempty(arg.texts)
		hlegend = legend(arg.texts);
		set(hlegend, 'Location', 'southoutside', 'FontSize', 8);
	end
	drawnow;
	hold off;
end

end

function time_series_clicker_test()
% warning('off', 'MATLAB:toeplitz:DiagonalConflict');
tmp = toeplitz(-2:4, -2:-2:-8);
img1 = cat(3, tmp, -1*tmp, 2*tmp, -2*tmp);
imgs = cat(4, img1, flipdim(img1, 1), flipdim(img1, 2), -img1);
texts = {'im1', 'im2', 'im3', 'im4'};
time_series_clicker(imgs, 'texts', texts);

end