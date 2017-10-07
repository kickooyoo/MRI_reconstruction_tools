function CAoA = concat_to_CAoA(CAoA, ndcs, stuff)
% concatenate to cell array of arrays (e.g, N x M cell array of inner arrays of variable length

contents = CAoA(ndcs);
aug_contents = cat(1, contents, stuff);

CAoA(ndcs) = aug_contents;
