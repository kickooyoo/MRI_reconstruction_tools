function reduced_fft = retrospective_undersample(mapped_im,SP)

% retrospectively undersamples simulated data
% inputs:
%   SP:         nx x ny logical, sampling pattern
%   mapped_im:  nx x ny x nc complex, for simulated data ???
%
% output:
%   reduced_fft: ns x nc

assert(size(SP,1)==size(mapped_im,1),'sizes of SP and mapped_im do not match');
assert(size(SP,2)==size(mapped_im,2),'sizes of SP and mapped_im do not match');

nc = size(mapped_im,3);
ns = sum(SP(:)); % number of samples kept (out of nx x ny)

reduced_fft = zeros(ns,nc);

for coil_ndx = 1:nc
    im_fft = fft2(mapped_im(:,:,coil_ndx));
    im_fft_vect = im_fft(:);
    SPv = SP(:);
    ndxSP = SPv'.*[1:length(SPv)];
    ndxSP = ndxSP(ndxSP~=0); % holds indices of masked pixels
    reduced_fft(:,coil_ndx) = im_fft_vect(ndxSP);
end
