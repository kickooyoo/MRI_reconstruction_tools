function n = get_ndcs(N)
%function n = get_ndcs(N)
% do integer indexing with guaranteed center at zero

if mod(N, 2) == 0
        n = [- N/2 : N/2 - 1].';
else
        n = [- (N - 1)/2 : (N - 1)/2].';
end

end

