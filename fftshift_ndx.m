function out = fftshift_ndx(in, N)
% gives fftshifted index back
%
% to do: N is odd?
% Mai Le

if in <= N/2
	out = N/2 + in;
else 
	out = in - N/2;
end

