function [phantom, motion] = toy_dynamic_phantom(varargin)
%function [phantom, motion] = toy_dynamic_phantom(varargin)
%

arg.Nx = 20;
arg.Ny = 22;
arg.Nz = 24;
arg.Nt = 30;
arg.radii = [5, 8]; 
arg.values = repmat(arg.Nt, 2); 
arg.resp_freq = 1/arg.Nt; % in samples
arg.card_freq = 5/arg.Nt; % in samples
arg.centers = [arg.Nx/2 arg.Ny/2 arg.Nz/2; arg.Nx/2 arg.Ny/4 arg.Nz + mean(arg.radii(:,2))/2];
arg.motion = []; % [Nt, numspheres, 3]
arg.cardiac = true;
arg = vararg_pair(arg, varargin);

num_spheres = size(arg.radii, 2);
assert(all(size(arg.centers) == [num_spheres 3]), 'invalid dimension for arg.centers');
assert(all(size(arg.values) == [arg.Nt num_spheres]), 'bad dims for values');

if isempty(arg.motion)
	motion = col(sin((1:arg.Nt)*2*pi*arg.resp_freq));
	arg.motion(:,:,3) = col(arg.radii/2)*motion;
else
	motion = col(mean(arg.motion(:,:,3),2));
end

if (size(arg.radii, 1) == 1) && (arg.Nt ~= 1)
        arg.radii = repmat(arg.radii, [arg.Nt 1]);
end
if arg.cardiac
        cmotion = col(sin((1:arg.Nt)*2*pi*arg.card_freq));
        arg.radii(:,1) = arg.radii(:,1).*cmotion;
end

if isinf(arg.resp_freq)
	arg.motion = zeros(arg.Nt, num_spheres, 3);
end

[xx, yy, zz] = ndgrid(1:arg.Nx, 1:arg.Ny, 1:arg.Nz);

for tt = 1:arg.Nt
	phantom(:,:,:,tt) = zeros(arg.Nx, arg.Ny, arg.Nz);
	for ss = 1:num_spheres
		curr_center = col(arg.centers(ss, :)) - col(arg.motion(tt, ss, :));
                curr_dist = sqrt((xx - curr_center(1)).^2 + (yy - curr_center(2)).^2 + (zz - curr_center(3)).^2);
		curr_sphere = arg.values(tt, ss)*(curr_dist < arg.radii(tt, ss));
		phantom(:,:,:,tt) = phantom(:,:,:,tt) + curr_sphere;
	end
end





