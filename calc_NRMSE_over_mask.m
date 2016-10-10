function NRMSE = calc_NRMSE_over_mask(x, xtrue, varargin)
% function NRMSE = calc_NRMSE_over_mask(x, xtrue, varargin)
% 
% calculates normalized root mean squared error between x and xtrue over mask

if nargin == 3
	mask = varargin{1};
else
	mask = true(size(x));
end

if numel(x) ~= numel(xtrue) || numel(x) ~= numel(mask)
	display('Error in calc_NRMSE_over_mask: sizes of inputs inconsistent.');
	keyboard;
end

xtrue = double(xtrue);
x = double(x);

%diffs = (abs(x(:))-abs(xtrue(:))).^2;
%diffs = abs(x(:)-xtrue(:)).^2;
%MSE = mean(diffs(logical(mask(:))));
MSE = norm(abs(x(:) - xtrue(:)), 2)^2;
NRMSE = sqrt(MSE)/norm(xtrue(:), 2);

