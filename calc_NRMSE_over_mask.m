function NRMSE = calc_NRMSE_over_mask(x,true,mask)
% function NRMSE = calc_NRMSE_over_mask(x,true,mask)
diffs = (abs(x)-abs(true)).^2;

MSE = mean(diffs(mask(:)));

NRMSE = sqrt(MSE)/max(abs(true(:)));
