function ph = gen_sim_phase(dims, phase_type)
% function ph = gen_sim_phase(dims, phase_type)
% input:
%       dims [Nx Ny (Nz)]
%       phase_type: 'ramp' or 'blob'
%
% to do: allow specialize constants

Nx = dims(1);
Ny = dims(2);
if length(dims == 3)
        Nz = dims(3);
end

switch phase_type
case 'ramp'
	[xx, yy, zz] = ndgrid(1:Nx, 1:Ny, 1:Nz);
	ph = 0.2*xx+0.5*yy-0.1*zz;
	ph = ph/max(col(ph))*0.9*pi;
case 'blob'
	tmp = mri_sensemap_sim_3D('nx', Nx, 'ny', Ny, 'nz', Nz, 'rcoil', Nx/2, 'ncoil', 8);
	tmp2 = abs(tmp);
	ph = 0.5*tmp2(:,:,:,1) - 0.7*tmp2(:,:,:,2) + 0.8*tmp2(:,:,:,3) - 0.2*tmp2(:,:,:,4) + ...
		0.5*tmp2(:,:,:,5) - 0.7*tmp2(:,:,:,6) + 0.8*tmp2(:,:,:,7) - 0.2*tmp2(:,:,:,8);
	ph = ph/max(col(ph))*0.9*pi;
otherwise
end