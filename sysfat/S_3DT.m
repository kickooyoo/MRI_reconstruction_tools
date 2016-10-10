function S = S_3DT(sense_maps, Nt, varargin)
%function S = S_3DT(sense_maps, Nt, varargin)
% generalized to dynamic sense_maps
%if nargin < 1, help(mfilename), error(mfilename), end
%if nargin == 1 && streq(varargin{1}, 'test'), Gdft_test, return, end
%
% object f should be Nx x Ny x Nz x Nt
assert(ismember(ndims_ns(sense_maps), [3 4]), 'invalid dimensions for sense_maps');
% if (ndims_ns(sense_maps) == 4) && (size(sense_maps,4) ~= Nt)
% 	% duplicate sense_maps for intermediate frames (rect basis)
% 	h = Nt/size(sense_maps,4);
% 	assert(mod(h,1) == 0, 'non integer ratio between Nt and sense_map frames');
% 	exp_sense_maps = permute(repmat(sense_maps, [1 1 1 1 h]), [1 2 3 5 4]);
% 	sense_maps = reshape(exp_sense_maps, size(sense_maps,1), size(sense_maps,2),...
% 		size(sense_maps,3), Nt);
% end
% 

arg.Nx = size(sense_maps, 1);
arg.Ny = size(sense_maps, 2);
arg.Nz = size(sense_maps, 3);
arg.Nr = arg.Nx * arg.Ny * arg.Nz;
arg.Nt = Nt;
arg.Nc = size(sense_maps, 4);
arg.mask = true(arg.Nx, arg.Ny, arg.Nz, arg.Nt);
arg = vararg_pair(arg, varargin);
arg.smaps = sense_maps;

if arg.Nc == 1
	S = Gdiag(ones(sum(col(arg.mask))*arg.Nt,1), 'mask', repmat(arg.mask, [1 1 arg.Nt]));
else
	if arg.Nz > 1
		if arg.Nt > 1
			S = fatrix2('idim', ([arg.Nx arg.Ny arg.Nz arg.Nt]) ,'arg',arg,'odim', ...
				([arg.Nx arg.Ny arg.Nz arg.Nt arg.Nc]), 'forw', @S_3DT_forw, ...
				'back', @S_3DT_back, 'imask', arg.mask);
		else
			S = fatrix2('idim', ([arg.Nx arg.Ny arg.Nz]) ,'arg',arg,'odim', ...
				([arg.Nx arg.Ny arg.Nz 1 arg.Nc]), 'forw', @S_3DT_forw, ...
				'back', @S_3DT_back, 'imask', arg.mask);

		end
	else
		if arg.Nt > 1
			S = fatrix2('idim', ([arg.Nx arg.Ny 1 arg.Nt]) ,'arg',arg,'odim', ...
				([arg.Nx arg.Ny 1 arg.Nt arg.Nc]), 'forw', @S_3DT_forw, ...
				'back', @S_3DT_back, 'imask', arg.mask);
		else
			S = fatrix2('idim', ([arg.Nx arg.Ny]) ,'arg',arg,'odim', ...
				([arg.Nx arg.Ny 1 1 arg.Nc]), 'forw', @S_3DT_forw, ...
				'back', @S_3DT_back, 'imask', arg.mask);
		end
	end
end

end

% y = G * x
function s = S_3DT_forw(arg,x)

%s = zeros(arg.Nx, arg.Ny, arg.Nz, arg.Nt, arg.Nc);
x_rep = repmat(x,[1 1 1 1 arg.Nc]);
smaps_perm = permute(arg.smaps, [1 2 3 5 4]);
S_rep = repmat(smaps_perm,[1 1 1 arg.Nt 1]);
s = x_rep.*S_rep;

end

% x = G' * y
function x = S_3DT_back(arg,s)

%x = zeros(arg.Nx, arg.Ny, arg.Nz, arg.Nt);

smaps_perm = permute(conj(arg.smaps), [1 2 3 5 4]);
S_rep = repmat(smaps_perm,[1 1 1 arg.Nt 1]);
prod_rep = S_rep.*s;
x = sum(prod_rep, 5);

end
