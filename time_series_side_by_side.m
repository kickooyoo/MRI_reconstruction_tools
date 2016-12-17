function time_series_side_by_side(time_series, varargin)
% function time_series_side_by_side(time_series, varargin)
%
% varargin:
%       yoffset
%       t
%	labels
% 	same_amp

if iscell(time_series)
        maxabs = @(t) max(abs(t));
        arg.yoffset = 2*max(cell2mat(cellfun(maxabs, time_series, 'UniformOutput', 0)));
else
        arg.yoffset = 2*max(abs(col(time_series)));
end
arg.t = [];
arg.labels = {};
arg.same_amp = false;
arg.draw_zero = false;
arg.FontSize = 10;
arg = vararg_pair(arg, varargin);

<<<<<<< Updated upstream
[Ntime, Nseries] = size(time_series);
if arg.same_amp
	amp = max(abs(col(time_series)));
	for ii = 1:Nseries 
		curr_max = max(abs(col(time_series(:,ii))));
		time_series(:,ii) = time_series(:,ii)*amp/curr_max;	
		if ~isempty(arg.labels)
			arg.labels{ii} = [arg.labels{ii} sprintf(', scale: %1.1e', curr_max)];
		end
	end
end

offsets = ones(Ntime, 1) * arg.yoffset*col(0:-1:-(Nseries - 1))';
x_offset = max(cellfun('length', arg.labels)) * arg.FontSize/3;
if isempty(arg.t)
        plot(offsets + time_series);
	hold on; plot(offsets, 'k--');
	text_x = Ntime - x_offset;
=======
if iscell(time_series) 
        figure; hold on;
        if iscell(arg.t)
                
                % to do: check dims of t
                
                Nentries = length(time_series);
                counter = 0;
                for ii = 1:Nentries
                        if ~isempty(time_series{ii})
                                if length(arg.t{ii}) ~= length(time_series{ii})
                                        keyboard
                                end
                                plot(arg.t{ii}, time_series{ii} + arg.yoffset*counter);
                                plot(arg.t{ii}, arg.yoffset*counter, 'k--');
                                counter = counter + 1;
                        end
                end
        
        elseif ~isempty(arg.t)
                % try to use vector t for all
                keyboard;
        else
                keyboard;
        end
>>>>>>> Stashed changes
else
        
        
        [Ntime, Nseries] = size(time_series);
        if arg.same_amp
                amp = max(abs(col(time_series)));
                for ii = 1:Nseries
                        curr_max = max(abs(col(time_series(:,ii))));
                        time_series(:,ii) = time_series(:,ii)*amp/curr_max;
                        arg.labels{ii} = [arg.labels{ii} sprintf(', scale: %1.1e', curr_max)];
                end
        end
        
        
        offsets = ones(Ntime, 1) * arg.yoffset*col(0:-1:-(Nseries - 1))';
        x_offset = max(cellfun('length', arg.labels)) * arg.FontSize/3;
        if isempty(arg.t)
                plot(offsets + time_series);
                hold on; plot(offsets, 'k--');
                text_x = Ntime - x_offset;
        else
                % to do: check dims t
                plot(arg.t, offsets + time_series)
                hold on; plot(arg.t, offsets, 'k--');
                text_x = max(arg.t) - x_offset;
        end

        
        
end

if ~isempty(arg.labels)
	text(double(text_x*ones(1,Nseries)), double(offsets(1,:) - arg.yoffset/2), arg.labels, 'FontSize', arg.FontSize, 'Interpreter', 'none');
end
