function true_image = make_sim_image(nx,ny)
%function true_image = make_sim_image(nx,ny)
%
% creates simulated static image with planar varying phase
% from brainweb, used in dynamic/static_SENSE_splitting project
% crops to nx x ny

assert(nx <= 258,'nx larger than maximum size of 258');
assert(ny <= 258,'ny larger than maximum size of 258');

fsep = filesep;
img_dir = [path_find_dir('mri') fsep '..' fsep 'data' fsep 'mri'];
img_path = [img_dir fsep 'brainweb_t1.jpg']; % 258x258
%img_dir = path_find_dir('tridiag_vega');
img_path = ['sag_brainweb_t1.jpg'];
true_image = double(imread(img_path));
dims = size(true_image);
[xx,yy] = ndgrid(-dims(1)/2:dims(1)/2-1,-dims(2)/2:dims(2)/2-1);
ph_max = pi/5;
ph = ph_max*xx/(dims(1)/2)+ph_max*yy/(dims(2)/2);
true_image = true_image.*exp(1i*ph);
extra_x = dims(1)-nx;
extra_y = dims(2)-ny;
if (mod(extra_x,2) == 0)
    true_image = true_image(1+extra_x/2:end-extra_x/2,:);
else
    true_image = true_image(1+(extra_x+1)/2:end-(extra_x-1)/2,:);
end
if (mod(extra_x,2) == 0)
    true_image = true_image(:,1+extra_y/2:end-extra_y/2);
else
    true_image = true_image(:,1+(extra_y+1)/2:end-(extra_y-1)/2);
end
