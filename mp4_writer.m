function mp4_writer(x, filename, varargin)
% function mp4_writer(x, filename, varargin)
% 
% inputs:
%       x [3D dyn image, 4D RGB dyn image, or cell of frame structs]
% 
% varargin: rate (fps)
%			magnify, 100 = regular
%			texts
%                       aspect
%
%       
% wrapper for Matlab's VideoWriter
arg.NI = 1; % number of images side by side
arg.rate = 1; % fps
arg.rgb = 0;
arg.magnify = 'fit'; % or can be positive number, 100 = regular size
arg.texts = {};
arg.text_x = 10*ones(1,arg.NI);
arg.text_y = 10 + size(x,2)/arg.NI*(0:arg.NI-1);%min(size(x,1), size(x,2))*(0:arg.NI-1);
arg.profile = 'MPEG-4';
arg.aspect = ones(1,3);
arg.debug = false;
if size(x, 1) < size(x, 2)
        tmp = arg.text_y;
        arg.text_y = arg.text_x;
        arg.text_x = tmp;
end
arg = vararg_pair(arg, varargin);

if strcmp(filename, 'tmp')
        if strcmp(arg.profile, 'MPEG-4')
                filename = '~/Downloads/tmp.mp4';
        elseif ~isempty(strfind(lower(arg.profile), 'avi'))
                filename = '~/Downloads/tmp.avi';
        else
                display('unknown file type');
                return;
        end
elseif length(filename) < 4 || ~strcmp(filename(end-3:end), '.mp4')
        if strcmp(arg.profile, 'MPEG-4')
                filename = [filename '.mp4'];
        elseif ~isempty(strfind(lower(arg.profile), 'avi'))
                filename = [filename '.avi'];
        else
                display('unknown file type');
                return;
        end
end


try
        writerObj = VideoWriter(filename, arg.profile);
catch
        display(sprintf('unable to write mp4 %s, this machine probably does not have the specified VideoWriter profile', filename))
        return;
end
writerObj.FrameRate = arg.rate;
open(writerObj);

if iscell(x)
        Nf = length(x);
else
        % normalize to satisfy writeVideo
        x = x./max(abs(col(x)));
        if arg.rgb
                Nf = size(x,4);
        else
                Nf = size(x,3);
        end
end
fhandle = figure;
for ii = 1:Nf
        if iscell(x)
                writeVideo(writerObj, x{ii})
        else
                if arg.rgb
                        pic = x(:,:,:,ii);
                else
                        pic = x(:,:,ii);
                end
                imshow(squeeze(pic), 'InitialMagnification', arg.magnify);
                daspect(arg.aspect);
                if ~isempty(arg.texts)
                        for jj = 1:length(arg.texts)
                                text(arg.text_x(jj), arg.text_y(jj), arg.texts{jj}, ...
                                        'color', [1 1 1]);
                        end
                end
                frame = getframe(fhandle);
                if arg.debug
                        cmap = colormap;
                        frame = im2frame(squeeze(pic)*255+1, cmap);
                end
                writeVideo(writerObj, frame);
        end
end
close(writerObj);
