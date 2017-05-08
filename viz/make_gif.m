function make_gif(x, filename, varargin)
%function make_gif(x, filename, varargin)

arg.png = true; % save series of png files for use with latex animate package
if isreal(x)
	arg.clim = [min(x(:)) max(x(:))];
else
	arg.clim = [min(abs(x(:))) max(abs(x(:)))];
end
arg.aspect = [1 1 1];
arg = vararg_pair(arg, varargin);

x = squeeze(x);

if ~isreal(x)
	x = abs(x);
	display('warning: taking absolute value of complex x');
end

Nf = size(x, 3);

if arg.png
	figure;
	for ff = 1:Nf
		im(x(:,:,ff,:), arg.clim);
		title('');
		colorbar;
		daspect(arg.aspect)
		drawnow;
		display(sprintf('printing png frame %d/%d', ff, Nf))
		print(sprintf('%s-%d',filename, ff), '-dpng');
	end
elseif ~exist([filename '.gif'])
	figure;
	for ff = 1:Nf
		im(x(:,:,ff), arg.clim);
		title('');
		colorbar;
		daspect(arg.aspect)
		drawnow;
		frame = getframe();
		img = frame2im(frame);
		[imind, cm] = rgb2ind(img,32);
		if ff == 1;
			imwrite(imind,cm, [filename '.gif'], 'gif', 'Loopcount',inf);
		else
			imwrite(imind,cm, [filename '.gif'],'gif','WriteMode','append');
		end
	end
else
	display(sprintf('filename %s.gif already exists!', filename));
end
end
