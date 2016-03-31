function [sense_maps, bodycoil_mask] = est_S_reg(coil_images, varargin)
%function [sense_maps, bodycoil_mask] = est_S_reg(coil_images, varargin)
% estimate sense maps from coil images
% basically wrapper for Michael's method
%
% inputs:
%	coil_images [Nx Ny Nc]
% 
% varargin:
%	figs_on
% 	thresh
%	l2b
%	dilate
% 
% parameters tuned for Siemes data 12/19
% 2D only :(
Nc = size(coil_images, 3);
arg.figs_on = 0;
arg.thresh = 0.2;% 0.45; %0.2 good for pat2
arg.l2b = 2;
arg.dilate = 3;% 5;% 3 good for pat2
arg.covmat = [];
arg.display_thresh = 1e-4;
arg.check_bodycoil_mask = false; % for quick pat2
arg.bodycoil = [];
arg.bodycoil_mask = [];
arg = vararg_pair(arg, varargin);

if isempty(arg.bodycoil)
	SoS = sos_combine(permute(coil_images, [1 2 4 3]), arg.covmat, []);
	bodycoil_sim = SoS;%.*exp(1i*angle(coil_images(:,:,9)));
	if 0
		SoS = fftshift(fft2(SoS));
		[nx ny] = size(SoS);
		figure; im(sOs);
		sOs_filt = sOs;
		window_width = 1;
		sOs_filt(nx/2-window_width:nx/2+window_width,ny/2-window_width:ny/2+window_width) = 0;
		SoS_filt = ifft2(ifftshift(sOs_filt));
		figure; im(SoS_filt);
	end
else
        SoS = arg.bodycoil;
         bodycoil_sim = SoS;
end
	% phase of coils was noisy, don't bother adding it
if isempty(arg.bodycoil_mask)
	arg.bodycoil_mask = adaptivethreshold(abs(SoS), 150, arg.thresh) & (abs(SoS) > arg.thresh/2*max(col(abs(SoS))));
	arg.bodycoil_mask = imerode(arg.bodycoil_mask, strel('disk', round(0.3*arg.dilate))); % get rid of extraneous pixels outside body
	arg.bodycoil_mask = bwconvhull(arg.bodycoil_mask);
	arg.bodycoil_mask = imerode(arg.bodycoil_mask, strel('disk', round(1.5*arg.dilate))); % further erode
	if sum(arg.bodycoil_mask) == 0
		display('empty body coil mask!');
		keyboard;
	end
end
bodycoil_mask = arg.bodycoil_mask;
bodycoil_sim = bodycoil_sim.*bodycoil_mask;
if arg.check_bodycoil_mask
        figure; im(bodycoil_sim)
        display('check masks, bodycoil_sim = redo_bodycoil_mask(SoS, arg)?');
        keyboard;
end

if ~all(size(bodycoil_sim) == size(coil_images(:,:,1)))
	display('sizes of bodycoil_sim and coil_images do not match');
	keyboard;
end
[sense_maps, sinit] = mri_sensemap_denoise(coil_images, 'bodycoil', bodycoil_sim, ...
	'chol', 1, 'niter', 1, 'l2b', arg.l2b);
if arg.figs_on
	% build mask for display
	mask = abs(bodycoil_sim)> arg.display_thresh;
	convex_mask = bwconvhull(mask); % only in 2014a, use [v d] = version;
	convex_mask = imdilate(convex_mask, strel('disk', arg.dilate));
	figure; im(sense_maps.*repmat(convex_mask, [1 1 Nc]));
end

end

function bodycoil_sim = redo_bodycoil_mask(SoS, arg)

bodycoil_sim = SoS;%.*exp(1i*angle(coil_images(:,:,9)));
bodycoil_mask = adaptivethreshold(SoS, 150, arg.thresh) & (SoS > arg.thresh/2*max(col(SoS)));
bodycoil_mask = imerode(bodycoil_mask, strel('disk', round(0.3*arg.dilate))); % get rid of extraneous pixels outside body
bodycoil_mask = bwconvhull(bodycoil_mask);
bodycoil_mask = imerode(bodycoil_mask, strel('disk', round(1.5*arg.dilate))); % further erode
bodycoil_sim = bodycoil_sim.*bodycoil_mask;
%figure; im(bodycoil_sim)

end
