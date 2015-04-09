function make_gif(x, filename)
%function make_gif(x, filename)

if ~isreal(x)
	x = abs(x);
	display('warning: taking absolute value of complex x');
end

if ~exist([filename '.gif'])
	figure;
	max_val = max(x);
	min_val = min(x);
	for ff = 1:size(x,3)
		im(x(:,:,ff), [min_val max_val]);
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
