function out = shift_cell_index(in, shift)
%function out = shift_cell_index(in, shift)
%
% out{ii} = in{ii + shift};

if mod(shift, 1) ~= 0
	display('invalid shift, must be integer')
	keyboard;
end

N = length(in);

for ii = 1:N
	if ii + shift > N
		curr_ii = mod(ii + shift, N);
		out{ii} = in{curr_ii};
	elseif ii + shift < 1
		curr_ii = mod(ii + shift, N);
		if curr_ii == 0
			curr_ii = N;
		end
		out{ii} = in{curr_ii};
	else
		out{ii} = in{ii + shift};
	end
end
