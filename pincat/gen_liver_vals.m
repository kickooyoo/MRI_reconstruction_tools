addpath('~/Documents/mai_code/spline_basis/methods/sim/dce/')

MR_params.r1 = 4.5; % for Gd at 3T, units mMol^-1 sec^-1
MR_params.flip_angle = 10/360*2*pi;
MR_params.rho0 = 1; % !!
% use true value for now, later need to use earlier T1 mapping scan to estimate this.
MR_params.TE = 0; % ignore for now
MR_params.T2star = 1; % ignore for now
MR_params.TR = 4.6 / 1000; % sec

MR_params.T1_read = 809 / 1000; % https://www.ncbi.nlm.nih.gov/pubmed/14990831
T1 = 809 / 1000;
T2 = 34 / 1000;

E1 = exp(-MR_params.TR./ T1);
% rho hat
true_image = MR_params.rho0.*sin(MR_params.flip_angle) .* div0(1-E1, 1 - E1*cos(MR_params.flip_angle)) .* ...
        exp(div0(-MR_params.TE, T2));

Nf = 10;
[C_aorta, C_t, C_p, duration, dt] = gen_AIF(Nf, 0.6, 2);
[dyn, MR_params] = add_contrast_agent(repmat(true_image, [1 1 Nf]), true, C_t, MR_params);

dyn_sc = dyn./true_image;

liver_base = 120;

d = round(liver_base * dyn_sc);
plot(d)

%% write to file for use by gen_params.py for pincat
write_fname = sprintf('liver_vals_%dNf.txt', Nf);
fileID = fopen(write_fname, 'w');
fprintf(fileID, '%d ', d);
fclose(fileID);

