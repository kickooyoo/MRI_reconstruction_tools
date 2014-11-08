function [img, full_kspace] = PF_3D(partial_kspace, full_dims, varargin)
%function [img, full_kspace] = PF_3D(partial_kspace, full_dims, varargin)
% 2D or 3D version of Xu and Haacke's 2001 "Partial Fourier Imaging in
% Multi-Dimensions: A Means to Save a Full Factor of Two in Time", JMRI 
% (xu:01:pfi)
%
% inputs:
% 	partial_kspace: [pf_Nx pf_Ny pf_Nz] 
%		assumes is taken from origin-centered k-space
% 	full_dims [Nx Ny Nz] : dims of full_kspace
%
% varargin:
%	PF_location: 
%		3x1 logical, placement of partial_kspace in full_kspace
%		0: 1:N
%		1: end-N+1:end
%		default: [1 0 1], for dynamic phantom data 2014/05/08
% 	filter_on: default = 1
%	niters: default = 5
%	window_step3 
%		width of transition band for apodization in low freq estimation
%		default = 8
%	window_step8 
%		width of transition band f);or apodization in iterative step 
%		default = 3
% 	show_conv_kern
%		to visualize how this method fills in missing kspace
%		default 0
%
% outputs:
%	img: [Nx Ny Nz] 
% 	kspace: [Nx Ny Nz]
%
% 2014-06-17 Mai Le
%
% notes: 
% 	- "filter" is more of an apodization, reduces discontinuities in k-space
% 	- can only handle even sizes for now
%

arg.PF_location = [1 0 1];
arg.filter_on = 1;
arg.niters = 5;
arg.window_step3 = 8;
arg.window_step8 = 3;
arg.show_conv_kern = 0;
arg = vararg_pair(arg, varargin);

% detect 2D or 3D
assert(all(mod(full_dims,2) == 0), 'err: odd dims for full kspace');
assert(ndims(partial_kspace) == length(full_dims),'full_dims does not match dims of partial_kspace');
if length(full_dims) == 2
	threeD = 0;
elseif length(full_dims) == 3
	threeD = 1;
else
	error('invalid number of dimensions, must be 2D or 3D');
end

% determine sizes of full, partial, and low-frequency kspace
Nx = full_dims(1);
Ny = full_dims(2);
[pf_Nx, pf_Ny, pf_Nz] = size(partial_kspace);
lf_Nx = (Nx/2 - (Nx - pf_Nx))*2;
lf_Ny = (Ny/2 - (Ny - pf_Ny))*2;
lf_xndx = (Nx - pf_Nx) + 1:pf_Nx;
lf_yndx = (Ny - pf_Ny) + 1:pf_Ny;
if threeD
	Nz = full_dims(3);
	lf_Nz = (Nz/2 - (Nz - pf_Nz))*2;
	lf_zndx = (Nz - pf_Nz) + 1:pf_Nz;
else
	Nz = 1;
	lf_Nz = 1;
	lf_zndx = 1;
end

% assign indeces for partial Fourier
pf_xndx = (1:pf_Nx) + arg.PF_location(1)*(Nx-pf_Nx);
pf_yndx = (1:pf_Ny) + arg.PF_location(2)*(Ny-pf_Ny);
pf_zndx = (1:pf_Nz) + arg.PF_location(3)*(Nz-pf_Nz);

% ------- begin algo -----------

% step 1: zero fill
full_kspace = zeros(Nx, Ny, Nz);
full_kspace(pf_xndx, pf_yndx, pf_zndx) = partial_kspace;
est_img = ifftn(ifftshift(full_kspace));

% step 2 and 3: Hanning filter on central k-space
lf_kspace = zeros(Nx, Ny, Nz);
lf_kspace(lf_xndx, lf_yndx, lf_zndx) = full_kspace(lf_xndx, lf_yndx, lf_zndx);
% selectively filter only edges of central k-space?

% build filter and masks for filter areas
if arg.filter_on
	hanning_dims_step3 = [length(lf_xndx) length(lf_yndx) length(lf_zndx)]; 
	hanning_step3 = generate_hanning_window(arg.window_step3, hanning_dims_step3, threeD);
	lf_kspace(lf_xndx, lf_yndx, lf_zndx) = lf_kspace(lf_xndx, lf_yndx, lf_zndx).*hanning_step3;
end

% step 4: estimate phase from lp image
est_lf_img = ifftn(ifftshift(lf_kspace));
phi = angle(est_lf_img); 

% build filter for step 8 outside for loop
if arg.filter_on
	hanning_dims_step8 = [pf_Nx pf_Ny pf_Nz] + arg.window_step8;
	hanning_step8 = generate_hanning_window(arg.window_step8, hanning_dims_step8, threeD);
	if threeD
		hanning_step8_crop = hanning_step8((1:pf_Nx) + arg.window_step8*(1-arg.PF_location(1)), (1:pf_Ny) + arg.window_step8*(1-arg.PF_location(2)), (1:pf_Nz) + arg.window_step8*(1-arg.PF_location(3)));
	else
		hanning_step8_crop = hanning_step8((1:pf_Nx) + arg.window_step8*(1-arg.PF_location(1)), (1:pf_Ny) + arg.window_step8*(1-arg.PF_location(2)));
	end
	hanning_step8_full = zeros(Nx, Ny, Nz);
	hanning_step8_full(pf_xndx, pf_yndx, pf_zndx) = hanning_step8_crop;
else
	hanning_step8_full = abs(full_kspace) > 0; 
end

for iter = 1:arg.niters
	% step 5: 
	rho_1 = abs(est_img).*exp(i*phi);

	if arg.show_conv_kern
		phi_diff = phi - angle(est_img);
		subplot(1,2,1);
		im(log10(abs(fftshift(fftn(exp(i*(phi_diff)))))+1e-10))
		subplot(1,2,2);
		im(abs(fftshift(fftn(exp(i*(phi_diff)))))>0)
		drawnow;
		pause(0.5)
	end

	% step 6: 
	s_1 = fftshift(fftn(rho_1));
	
	% step 7 & 8: enforce data consistency but smoothly
	tmp = hanning_step8_full.*full_kspace + (1-hanning_step8_full).*s_1;
	
	% step 9:
	est_img = ifftn(ifftshift(tmp));
end

img = est_img;
full_kspace = tmp;

end

function hanning_window = generate_hanning_window(window_length, dims, threeD)
% for apodization, window length refers to transition band length
% dims is size of window
	if any(dims(1:end-1*(1-threeD)) <= window_length*2)
		display('ERR: transition band wider than provided dims for Hanning window');
		keyboard;
	end

	% add extra 3: 2 because end points are zero, 1 in the middle for expanding
	hann_len = (window_length+1)*2 + 1;
	hanning_1D = col(hann(hann_len));
	hanning_2D = hanning_1D*hanning_1D';
	if threeD
		hanning_3D = repmat(hanning_2D,[1 1 hann_len]).*repmat(permute(hanning_1D,[2 3 1]),[hann_len hann_len 1]);
		hanning_crop = hanning_3D(2:end-1,2:end-1,2:end-1)./max(col(hanning_3D)); 
	else
		hanning_crop = hanning_2D(2:end-1,2:end-1)./max(col(hanning_2D)); 
	end

	hanning_stretch = cat(1, hanning_crop(1:window_length,:,:), ...
		repmat(hanning_crop(window_length+1,:,:), [dims(1)-2*window_length 1 1]), ...
		hanning_crop(window_length+2:end,:,:));
	hanning_stretch_2D = cat(2, hanning_stretch(:,1:window_length,:), ...
		repmat(hanning_stretch(:,window_length+1,:), [1 dims(2)-2*window_length 1]), ...
		hanning_stretch(:,window_length+2:end,:));

	if threeD
		hanning_stretch_3D = cat(3, hanning_stretch_2D(:,:,1:window_length), ...
			repmat(hanning_stretch_2D(:,:,window_length+1), [1 1 dims(3)-2*window_length]), ...
			hanning_stretch_2D(:,:,window_length+2:end));
		hanning_window = hanning_stretch_3D;
	else
		hanning_window = hanning_stretch_2D;
	end
end
