function [avg_Hz, top_Hz, top_weights] = average_frequency(signal, sampling_period, varargin)
%function avg_Hz = average_frequency(signal, sampling_period)
% bit of a misnomer, really gives weighted average of top harmonics
% freqencies binned with nonzero tolerance
%
% inputs:
%	signal [1 N] (real-valued)
% 	sampling period (double) (seconds)
%
% vararg:
%	harmonic_tolerance 
%		percent error for gropuing together like frequencies and harmonics 
%		default: 0.25
% outputs:
%	avg_Hz (double) (seconds)
%		weighted average fundamental harmonic
%		avg_Hz = 1/sum(top_weights) * sum(top_Hz .* top_weights)
%	top_Hz	[1 M] (double)
%		top fundamental harmonics in descending order
%	top_weights [1 M] (double)
%		energy in frequency domain associated with each harmonic

arg.harmonic_tolerance = 0.25; % set o zero for no harmonic grouping at all
arg = vararg_pair(arg, varargin);

N = length(signal);
X = fft(signal);
freqs = (0:N - 1)/(sampling_period*N);

energy = abs(X).^2;

% look at first half of spectrum
energy = energy(1:ceil(N/2));
freqs = freqs(1:ceil(N/2));

% group together harmonics
[sort_weights, sort_order] = sort(energy,'descend');
sort_freqs = freqs(sort_order);

harmonics = [];
harmonic_weights = [];
for ii = 1:ceil(N/2)
	curr_freq = sort_freqs(ii);
	if isempty(harmonics)
		upper_harmonics = harmonics;
	else
		upper_harmonics = (1:4)'*harmonics;
	end
	match = abs(upper_harmonics - curr_freq)./curr_freq < arg.harmonic_tolerance;
	if sum(match) == 0
		harmonics = cat(2, harmonics, curr_freq);
		harmonic_weights = cat(2, harmonic_weights, sort_weights(ii));
	else
		if sum(match) > 1
			display('add weight to multiple? split it?');
			keyboard
		end
		bin = find(sum(match,1));
		harmonic_weights(bin) = harmonic_weights(bin) + sort_weights(ii);
	end
end

%figure; scatter(harmonics, harmonic_weights); title(sprintf('tol %.2f', arg.harmonic_tolerance))

avg_Hz = harmonics*harmonic_weights'/sum(harmonic_weights);

[top_weights, sort_harm_order] = sort(harmonic_weights, 'descend');
top_Hz = harmonics(sort_harm_order);



