function mp4_writer(x, filename, varargin)
% wrapper for Matlab's VideoWriter
arg.rate = 1; % fps
arg = vararg_pair(arg, varargin);

if ~strcmp(filename(end-3:end), '.mp4')
	filename = [filename '.mp4'];
end

% normalize to satisfy writeVideo
x = x./max(abs(col(x)));

writerObj = VideoWriter(filename, 'MPEG-4');
writerObj.FrameRate = arg.rate;
open(writerObj);
for ii = 1:size(x,3)
	im(x(:,:,ii));
	frame = getframe;
	writeVideo(writerObj, frame);
end
close(writerObj);