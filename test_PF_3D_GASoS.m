% script for testing PF_3D on simulated GA stack-of-stars data

if ~exist('phantom3d')
	GRASP_setup;
end
Nx = 128;
Ny = 100;
Nz = 60;
img = phantom3d('Modified Shepp-Logan', max([Nx Ny Nz]));
img = crop_outside(img, [Nx Ny Nz]); 

Nro = round(Nx*1);
Nspokes = 500;
Nslice = Nz;
Nslice_PF = round(0.8*Nslice);
grad_shift = 0; % 0.75?
k = get_GA_coords(Nro, Nspokes, Nslice, 'grad_shift', grad_shift);

Ns = numel(k);
Nt = 1;
Nc = 1;
F = F_NC_3DT(k, Ns, Nro, Nt, Nx, Ny, Nz, Nc);



