function output = resize3D(input, new_dims)
%function output = resize3D(input, new_dims)
% expect new_dims to be 3D
% preserve size of fourth dim
% leverages imresize, so uses bicubic interpolation in x-y then x-z

for ii = 1:size(input,3)
	for jj = 1:size(input,4)
		tmp(:,:,ii,jj) = imresize(input(:,:,ii,jj), new_dims(1:2));
	end
end
for ii = 1:new_dims(2)
	for jj = 1:size(input,4)
		output(:,ii,:,jj) = squeeze(imresize(squeeze(tmp(:,ii,:,jj)), new_dims([1 3])));
	end
end
