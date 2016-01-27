function NRMSE = calc_NRMSE_over_mask(x, xtrue, mask)
% function NRMSE = calc_NRMSE_over_mask(x, xtrue, mask)
% 
% calculates normalized root mean squared error between x and xtrue over mask

if numel(x) ~= numel(xtrue) || numel(x) ~= numel(mask)
	display('Error in calc_NRMSE_over_mask: sizes of inputs inconsistent.');
	keyboard;
end

diffs = (abs(x(:))-abs(xtrue(:))).^2;

MSE = mean(diffs(logical(mask(:))));

NRMSE = sqrt(MSE)/max(abs(xtrue(:)));
