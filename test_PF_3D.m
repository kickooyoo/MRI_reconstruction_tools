% test PF_3D.m with Shepp-Logan to examine behavior of xu:01:pfi

%% --- Parameters that control recon quality

percentage = 0.65; % partial k-space fraction, must be greater than 0.5

window_step3 = 12; % controls transition band width for apodization
window_step8 = 6;


%% --- Chose test case

Nz = 1;
% uncomment line below if want to test third dim
%Nz = 2*ceil((2*max(window_step3,window_step8)+1)/(2*(percentage - 0.5))/2); % even number just big enough to pass condition in PF_3D that requires apodization area of 16 

% Shepp-Logan
p = repmat(phantom('Modified Shepp-Logan'),[1 1 Nz]).*exp(-i*0.05); % because need even dims
[Nx Ny ~] = size(p);
k = fftshift(fftn(p));  

% Tough Case
%p = rand(256, 256, Nz).*(exp(i*0.05)*ones(256,256,Nz)); % not quite constant phase :(
%[Nx Ny ~] = size(p);
%k = fftshift(fftn(p));  

% Artificial k-space for visualizing effect of algo
%k = rand(size(p));
%p = ifftn(ifftshift(k));

%% --- Select partial k-space and do PF recon

partial_kspace = k(1:round(percentage*Nx),end-round(percentage*Ny)+1:end,1:round(percentage*Nz));
sampling = zeros(size(p));
sampling(1:round(percentage*Nx),end-round(percentage*Ny)+1:end,1:round(percentage*Nz)) = 1;
PF_location = [0 1 0];
if Nz > 1
	full_dims = [Nx Ny Nz];
else
	full_dims = [Nx Ny];
end
[img, full_kspace] = PF_3D(partial_kspace, full_dims, 'PF_location', PF_location, 'filter_on', 1, 'niters',20,'window_step3',window_step3,'window_step8',window_step8,'show_conv_kern',1);
diff_img = img-p;

figure; 
subplot(2,2,1); im(p)
subplot(2,2,2); im(img);
subplot(2,2,3); im(sampling);
subplot(2,2,4); im(diff_img); 

figure; 
klog = log10(abs(k)+1e-10);
subplot(1,2,1); im(klog);
subplot(1,2,2); im(log10(abs(full_kspace)+1e-10), [min(col(klog)) max(col(klog))]);

% conclusions: 
% - perfect recovery not possible 
% - missing corners of k-space filled in by convolution of fft of phase difference between current iterate and low frequency estimate

