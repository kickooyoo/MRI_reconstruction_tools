function dcf_tweak = tweak_dcf(dcf, freqs)
% function dcf_tweak = tweak_dcf(dcf, freqs)
% 
% heuristic tweaking of dcf 
dcf_tweak = max(0, dcf - 0.1e-6*exp(-abs(freqs)));

% figure; scatter(abs(freqs(1:5:end)), dcf(1:5:end));