function convexImage = convex_image(image)
% convexImage = function convex_image(image)
% 
% regionprops ConvexImage frustratingly returns cropped convex image
% this provides convexImage with same footprint as input image

imsize = size(image);
nd = ndims(image);
info = regionprops(image, 'ConvexImage', 'BoundingBox');
boxsize = info.BoundingBox(end/2+1:end);
boxstart = ceil(info.BoundingBox(1:end/2));
convexImage = false(imsize);
convexImage(boxstart(2) : boxstart(2) + boxsize(2) - 1, boxstart(1) : boxstart(1) + boxsize(1) - 1) = info.ConvexImage;

% figure; imagesc(image - convexImage);
