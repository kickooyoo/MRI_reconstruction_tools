function filt_out = Hamming_bandpass(in, freq1_Hz, freq2_Hz, sampling_period)
%function filt_out = Hamming_bandpass(in, freq1_Hz, freq2_Hz, sampling_period)
% inputs
%	sampling_period [seconds]
% for spincemaille:11:zip, window is 0.04 to 0.05 Hz
% in [Ntime Nsignals] for sampling multiple signals at once
% todo parallelize with filter2

order = 48;

[Ntime, Nsignals] = size(in);

% 2 to make it a fraction of the Nyquist rate
norm_freq1 = Hz_to_norm_freq(freq1_Hz, sampling_period)*2;
norm_freq2 = Hz_to_norm_freq(freq2_Hz, sampling_period)*2;
b = fir1(order, [norm_freq1 norm_freq2], 'bandpass');

for ii = 1:Nsignals
	curr_in = squeeze(in(:,ii));
	filt_out(:,ii) = conv(curr_in, b, 'same');
end

end

function norm_freq = Hz_to_norm_freq(Hz, sampling_period)

norm_freq = Hz*sampling_period; % 1/sec * sec/sample  = 1/sample

end
