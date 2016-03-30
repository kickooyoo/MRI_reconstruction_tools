function val = ind2sub2val(A, size, indeces)
% function val = ind2sub2val(A, size, indeces)
% extension of ind2sub, then gets value from matrix A belonging to sub
% 3D
[new_x, new_y, new_z] = ind2sub(size, indeces);
for ii = 1:length(indeces)
        val(ii) = A(new_x(ii), new_y(ii), new_z(ii));
end


