function val = ind2sub2val(A, size, indices)
% function val = ind2sub2val(A, size, indices)
% extension of ind2sub, then gets value from matrix A belonging to sub
% 3D
[new_x, new_y, new_z] = ind2sub(size, indices);
for ii = 1:length(indices)
        val(ii) = A(new_x(ii), new_y(ii), new_z(ii));
end


