function [phantom, motion] = toy_dynamic_phantom(varargin)
%function [phantom, motion] = toy_dynamic_phantom(varargin)
%

arg.Nx = 20;
arg.Ny = 22;
arg.Nz = 24;
arg.Nt = 30;
arg.radii = [5, 8];
arg.values = [1 5];
arg.centers = [arg.Nx/2 arg.Ny/2 arg.Nz/2; arg.Nx/2 arg.Ny/2 arg.Nz + arg.radii(2)/2];
arg.motion = []; % [num_spheres, Nt, 3]
arg = vararg_pair(arg, varargin);

num_spheres = length(arg.radii);
assert(all(size(arg.centers) == [num_spheres 3]), 'invalid dimension for arg.centers');

if isempty(arg.motion)
	motion = col(sin((1:arg.Nt)*pi*4/arg.Nt));
	arg.motion(:,:,3) = col(arg.radii/2)*permute(motion, [2 1 3]);
end

[xx, yy, zz] = ndgrid(1:arg.Nx, 1:arg.Ny, 1:arg.Nz);

for tt = 1:arg.Nt
	phantom(:,:,:,tt) = zeros(arg.Nx, arg.Ny, arg.Nz);
	for ss = 1:num_spheres
		curr_center = col(arg.centers(ss, :)) - col(arg.motion(ss, tt, :));
		curr_sphere = arg.values(ss)*(sqrt((xx - curr_center(1)).^2 + (yy - curr_center(2)).^2 + (zz - curr_center(3)).^2) < arg.radii(ss));
		phantom(:,:,:,tt) = phantom(:,:,:,tt) + curr_sphere;
	end
end





