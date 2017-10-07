function abc = make_alphabet(varargin)
% function abc = make_alphabet(varargin)
arg.range = 1:26;
arg.dims = [128 128];
arg.maxval = 1;
arg = vararg_pair(arg, varargin);

position = [floor(arg.dims(1)/10) floor(arg.dims(2)/3)];
FontSize = min(ceil(0.8*arg.dims));

alphabet = upper('abcdefghijklmnopqrstuvwxyz');
figure;
for ii = arg.range
        tmp = zeros(arg.dims);
        imshow(tmp);
        axis off; 
        text(position(1), arg.dims(2) - position(2), alphabet(ii), 'Color', 'white', 'FontSize', FontSize);
        frame = getframe(gca);
        frame = double(frame.cdata(:,:,1));
        abc(:,:,ii) = arg.maxval*frame/max(col(frame));
end
abc = abc(1:arg.dims(1), 1:arg.dims(2), :);
close