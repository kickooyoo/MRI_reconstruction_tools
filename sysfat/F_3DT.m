function F = F_3DT(Nx, Ny, Nz, Nt, Nresp, Nc, varargin)
%function F = F_3DT(Nx, Ny, Nz, Nt, Nresp, Nc,i varargin)
%
%if nargin < 1, help(mfilename), error(mfilename), end
%if nargin == 1 && streq(varargin{1}, 'test'), Gdft_test, return, end
%
% object f should be Nx x Ny x Nz x Nt x Nresp x Nc
% varargin: samp for subsampling (dynamic)
%
% obsolete! use F_5D instead
% to do: revert back to only 3D + time
display('currently not passing test_adjoint! use F_5D instead')

arg.Nx = Nx;
arg.Ny = Ny;
arg.Nz = Nz;
arg.Nt = Nt;
arg.Nresp = Nresp;
arg.Nc = Nc;

arg.Nr = arg.Nx*arg.Ny*arg.Nz;
arg.Ns = arg.Nr;
arg.samp = []; 
arg = vararg_pair(arg, varargin);

if ~isempty(arg.samp) && ( (size(arg.samp,1) ~= Nx) || (size(arg.samp,2) ~= Ny) || (size(arg.samp,3) ~= Nz))
        display(sprintf('samp size [%d %d %d] does not match given dims:[%d %d %d]', ...
                size(arg.samp,1), size(arg.samp,2), size(arg.samp,3), Nx, Ny, Nz));
        keyboard
end
if ~isempty(arg.samp) && ~islogical(arg.samp)
        display('samp must be logical');
        keyboard;
end
if ~isempty(arg.samp)
        arg.Ns = sum(col(arg.samp));
end

if isempty(arg.samp)
        F = fatrix2('idim', squeeze([arg.Nx arg.Ny arg.Nz arg.Nt arg.Nresp arg.Nc]), 'arg', arg, ...
                'odim', squeeze([arg.Nx arg.Ny arg.Nz arg.Nt arg.Nresp arg.Nc]), 'forw', ...
                @F_3DT_forw, 'back', @F_3DT_back);
else
        F = fatrix2('idim', squeeze([arg.Nx arg.Ny arg.Nz arg.Nt arg.Nresp arg.Nc]), 'arg', arg, ...
                'odim', [arg.Ns arg.Nt], 'forw', @F_3DT_forw, ...
                'back', @F_3DT_back);
        
end

end

% y = G * x
function S = F_3DT_forw(arg,s)

S = fft(fft(fft(s,[],1),[],2),[],3);
if ~isempty(arg.samp)
        S = S(repmat(arg.samp, [1 1 1 1 arg.Nc arg.Nresp]));
end

end

% x = G' * y
function s = F_3DT_back(arg,S)

if ~isempty(arg.samp)
        S = embed(col(S), repmat(arg.samp, [1 1 1 1 arg.Nc arg.Nresp]));
end
s = ifft(ifft(ifft(S,[],1),[],2),[],3)*arg.Nr;

end
