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
%	method (str) {'harmonics', 'autocorr'};
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
arg.method = 'autocorr';%'harmonics'; 
arg.harmonic_tolerance = 0.1; % set o zero for no harmonic grouping at all
arg.num_harmonics = 8;%4;
arg.freq_window = [0.1 1]; % use for autocorr
arg = vararg_pair(arg, varargin);

N = length(signal);
X = fft(signal);

switch arg.method
case 'harmonics'
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
			upper_harmonics = (1:arg.num_harmonics)'*harmonics;
		end
		bin_close = abs(upper_harmonics - curr_freq)./curr_freq;
		match = bin_close < arg.harmonic_tolerance;
		if sum(match) == 0
			harmonics = cat(2, harmonics, curr_freq);
			harmonic_weights = cat(2, harmonic_weights, sort_weights(ii));
		else
			best_bin = ceil(find(bin_close == min(col(bin_close)))/arg.num_harmonics);
			if length(best_bin) > 1
				%display('add weight to multiple? split it?');
				best_bin = best_bin(1);
			end
			%bin = find(sum(match,1));
			harmonic_weights(best_bin) = harmonic_weights(best_bin) + sort_weights(ii);
		end
	end

	%figure; scatter(harmonics, harmonic_weights); title(sprintf('tol %.2f', arg.harmonic_tolerance))

	avg_Hz = harmonics*harmonic_weights'/sum(harmonic_weights);

	[top_weights, sort_harm_order] = sort(harmonic_weights, 'descend');
	top_Hz = harmonics(sort_harm_order);

case 'autocorr'

% 	lags = 5:min(20, length(signal));
        lags = 1./(arg.freq_window * sampling_period);
        lags = max(round(lags(2)), 1) : min(round(lags(1)), length(signal));
	for ii = 1:length(lags)
		if lags(ii) < 0
			autoc(ii) = mean(abs(signal(1:end + lags(ii)) - signal(1-lags(ii):end)).^2);
		else
	 		autoc(ii) = mean(abs(signal(1+lags(ii):end)-signal(1:end -lags(ii))).^2);
		end
        end
        [autoc_peaks, autoc_peak_ndcs] = findpeaks(autoc, 'MinPeakDistance', lags(1));
	period = lags(autoc_peak_ndcs(abs(autoc_peaks - max(autoc)) < 1e-2)); % samples
	if length(period) ~= 1
                period = mode(period - [0 period(1:end-1)]);
        end
	if abs(period) <  min(lags), keyboard; end
	avg_Hz = 1/(period * sampling_period);% 1 / (sample/cycle * seconds/sample) = cycles/second = Hz
	top_weights = []; % to do: fill in runner up lags
	top_Hz = [];
% 	figure; plot(signal); hold on; plot(signal(1+period:end), 'r'); title(sprintf('period: %d samp, %1.1f Hz', period, avg_Hz));
% 	keyboard
otherwise 
	display(sprintf('unknown option %s', arg.method));
	keyboard
end


