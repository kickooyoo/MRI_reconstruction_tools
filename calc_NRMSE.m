function NRMSE = calc_NRMSE(x, xtrue, varargin)
% function NRMSE = calc_NRMSE(x, xtrue, varargin)
% 
% calculates normalized root mean squared error between x and xtrue over mask

if nargin == 1 && strcmp(x, 'test')
        test_calc_NRMSE();
        return
elseif nargin == 3
	mask = varargin{1};
else
	mask = true(size(x));
end

if numel(x) ~= numel(xtrue) || numel(x) ~= numel(mask)
	display('Error in calc_NRMSE: sizes of inputs inconsistent.');
	keyboard;
end
if nargin > 2
        if numel(x) ~= numel(mask)
                display('Error in calc_NRMSE: sizes of inputs inconsistent.');
                keyboard;
        end
        x = x(mask(:));
        xtrue = xtrue(mask(:));
end

xtrue = double(xtrue);
x = double(x);

MSE = norm(abs(x(:) - xtrue(:)), 2)^2;
NRMSE = sqrt(MSE)/norm(xtrue(:), 2);
end


function test_calc_NRMSE()
x = phantom();
xnoisy = x + randn(size(x));
mask = convex_image(x > 0);
figure; subplot(1,3,1); imagesc(x); title('xtrue'); 
subplot(1,3,2); imagesc(xnoisy); title('x')
subplot(1,3,3); imagesc(mask); title('mask')
NRMSE_mask = calc_NRMSE(xnoisy, x, mask)
NRMSE = calc_NRMSE(xnoisy, x)
end