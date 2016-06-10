% script for testing PF_3D on simulated GA stack-of-stars data

if ~exist('phantom3d')
	GRASP_setup;
end
Nx = 128;
Ny = 100;
Nz = 60;
img = phantom3d('Modified Shepp-Logan', max([Nx Ny Nz]));
img = crop_outside(img, [Nx Ny Nz]); 
phase_type = 'blob';
ph = gen_sim_phase([Nx Ny Nz], phase_type);

img = img.*exp(1i*ph);

Nro = round(Nx*1);
Nspokes = 500;
Nslice = Nz;
Nslice_PF = round(0.8*Nslice);
Nslice_nocollect = Nslice - Nslice_PF;
grad_shift = 0.75; % 0.75?
k = get_GA_coords(Nro, Nspokes, 1, 'grad_shift', grad_shift);

Ns = Nro*Nspokes;
Nt = 1;
Nc = 1;
F = F_NC_3DT(k, Ns, Nro, Nt, Nx, Ny, Nz, Nc);

rawdat = F*img;

rawdat = reshape(rawdat, Nro, Nspokes, Nz);
rawdat_PF = rawdat(:, :, end - Nslice_PF + 1: end);
dcf = repmat(radial_dcf(Nro), [1 Nspokes Nslice]);
Frawdat = F'*col(dcf.*rawdat);
Frawdat = reshape(Frawdat, Nx, Ny, Nz);

%% use 3D PF
[~, k_3D] = PF_3D(rawdat_PF, [Nro Nspokes Nslice], 'PF_location', [0 0 1], 'window_step3', 3);

img_3D = F'*col(dcf.*k_3D);
img_3D = reshape(img_3D, Nx, Ny, Nz);

%% use 3D PF multicoil
[k_3D_MC, ~] = PF_3D_MC(permute(rawdat_PF, [1 4 2 3]), 'full_dims', [Nro Nspokes Nslice]);

k_3D_MC = squeeze(k_3D_MC);
img_3D_MC = F'*col(dcf.*k_3D_MC);
img_3D_MC = reshape(img_3D_MC, Nx, Ny, Nz);

%% use 1D PF
% resulting in scaling! :(
%addpath('~/Documents/mai_code/partial_kspace');
overlap = Nslice_PF - Nslice/2;
for spoke_ndx = 1:Nspokes
	[~, ~, ~, k_1D(:, spoke_ndx, :)] = homodyne_recon(squeeze(rawdat_PF(:, spoke_ndx, :)), Nro, Nslice, overlap, 'direction', 2);
end
img_1D = F'*col(dcf.*k_1D);
img_1D = reshape(img_1D, Nx, Ny, Nz);

norm(col(img_3D - Frawdat))
norm(col(img_3D_MC - Frawdat))
norm(col(img_1D - Frawdat))
