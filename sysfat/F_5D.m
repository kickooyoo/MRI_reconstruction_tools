function F = F_5D(Nx, Ny, Nz, Nt, Nresp, Nc, varargin)
%function F = F_5D(Nx, Ny, Nz, Nt, Nresp, Nc, varargin)
%
%if nargin < 1, help(mfilename), error(mfilename), end
%if nargin == 1 && streq(varargin{1}, 'test'), Gdft_test, return, end
%
% object f should be Nx x Ny x Nt
% varargin: samp for subsampling (dynamic) [Nx Ny Nt Nresp] (assume same in Nz and Nc)
% to do: uneven numbre of samples for frames, resp states

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

if ~isempty(arg.samp) && ( (size(arg.samp,1) ~= Nx) || (size(arg.samp,2) ~= Ny) || (size(arg.samp,3) ~= Nz) )
        display(sprintf('samp size [%d %d %d] does not match given dims:[%d %d %d]', size(arg.samp,1), size(arg.samp,2), size(arg.samp,3), Nx, Ny, Nz));
        keyboard
end
if ~isempty(arg.samp) && ~islogical(arg.samp)
        display('samp must be logical');
        keyboard;
end
if ~isempty(arg.samp)
        arg.Ns = sum(col(arg.samp));
end

if (arg.Nc > 1)
        if isempty(arg.samp)
                F = fatrix2('idim', [arg.Nx arg.Ny arg.Nz arg.Nt arg.Nresp arg.Nc], 'arg', ...
                        arg, 'odim', [arg.Nx arg.Ny arg.Nz arg.Nt arg.Nresp arg.Nc], 'forw', ...
                        @F_5D_forw, 'back', @F_5D_back);
        else
                F = fatrix2('idim', [arg.Nx arg.Ny arg.Nz arg.Nt arg.Nresp arg.Nc], 'arg', ...
                        arg, 'odim', [arg.Ns arg.Nt arg.Nc], 'forw', ... % ??? odims
                        @F_5D_forw, 'back', @F_5D_back);
                
        end
else
        if isempty(arg.samp)
                F = fatrix2('idim', [arg.Nx arg.Ny arg.Nz arg.Nt arg.Nresp], 'arg', arg, ...
                        'odim', [arg.Nx arg.Ny arg.Nz arg.Nt arg.Nresp], 'forw', ...
                        @F_5D_forw, 'back', @F_5D_back);
        else
                F = fatrix2('idim', [arg.Nx arg.Ny arg.Nz arg.Nt arg.Nresp], 'arg', arg, ...
                        'odim', [arg.Ns arg.Nt], 'forw', @F_5D_forw, ... % ??? odims
                        'back', @F_5D_back);
                
        end
end

end

% y = G * x
function S = F_5D_forw(arg,s)

% % Fourier Transform
% S = zeros(arg.Nx,arg.Ny,arg.Nt,arg.Nc);
% for tbasis_ndx = 1:arg.Nt
% 	for coil_ndx = 1:arg.Nc
% 		S(:,:,tbasis_ndx,coil_ndx) = fft2(s(:,:,tbasis_ndx,coil_ndx)); 
% 	end
% end

S = fft(fft(fft(s,[],1),[],2),[],3);
if ~isempty(arg.samp)
        S = S(repmat(permute(arg.samp, [1 2 5 3 4 6]), [1 1 arg.Nz 1 1 arg.Nc]));
end

end

% x = G' * y
function s = F_5D_back(arg,S)

% s = zeros(arg.Nx,arg.Ny,arg.Nt,arg.Nc);
% for tbasis_ndx = 1:arg.Nt
%     for coil_ndx = 1:arg.Nc
% 	s(:,:,tbasis_ndx,coil_ndx) = arg.Nr*ifft2(S(:,:,tbasis_ndx,coil_ndx));
%     end
% end
if ~isempty(arg.samp)
        S = embed(col(S), repmat(permute(arg.samp, [1 2 5 3 4 6]), [1 1 arg.Nz 1 1 arg.Nc]));
end
s = ifft(ifft(ifft(S,[],1),[],2),[],3)*arg.Nr;

end
