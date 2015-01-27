function sense_maps = est_S_reg(coil_images, varargin)
% estimate sense maps from coil images
% basically wrapper for Michael's method
% 
% parameters tuned for Siemes data 12/19
% 2D only :(
arg.figs_on = 1;
arg.thresh = 0.15;
arg.l2b = 2;
arg.dilate = 10;
arg = vararg_pair(arg, varargin);
Nc = size(coil_images, 3);

SoS = sos_combine(permute(coil_images, [1 2 4 3]),[],[]);
bodycoil_sim = SoS;%.*exp(1i*angle(coil_images(:,:,9)));
% phase of coils was noisy, don't bother adding it
bodycoil_sim = bodycoil_sim.*(SoS > arg.thresh*max(col(SoS)));
[sense_maps, sinit] = mri_sensemap_denoise(coil_images, 'bodycoil', bodycoil_sim, ...
	'chol', 1, 'niter', 1, 'l2b', arg.l2b);
if arg.figs_on
	% build mask for display
	mask = abs(bodycoil_sim)>1e-4;
	convex_mask = bwconvhull(mask); % only in 2014a, use [v d] = version;
	convex_mask = imdilate(convex_mask, strel('disk', arg.dilate));
	figure; im(sense_maps.*repmat(convex_mask, [1 1 Nc]));
end

