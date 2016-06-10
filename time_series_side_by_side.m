function time_series_side_by_side(time_series, varargin)
% function time_series_side_by_side(time_series, varargin)
%
% varargin:
%       yoffset
%       t
arg.yoffset = 2*max(abs(col(time_series)));
arg.t = [];
arg.same_amp = false;
arg = vararg_pair(arg, varargin);

[Ntime, Nseries] = size(time_series);
if arg.same_amp
	amp = max(abs(col(time_series)));
	for ii = 1:Nseries 
		curr_max = max(abs(col(time_series(:,ii))));
		time_series(:,ii) = time_series(:,ii)*amp/curr_max;
	end
end

offsets = ones(Ntime, 1) * arg.yoffset*col(0:-1:-(Nseries - 1))';
if isempty(arg.t)
        plot(offsets + time_series)
else
        plot(arg.t, offsets + time_series)
end
