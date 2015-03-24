function make_gif(x, filename)
%function make_gif(x, filename)

if ~exist([filename '.gif'])
	figure; 
		    for ff = 1:size(x,3)
   			im(x(:,:,ff));
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
