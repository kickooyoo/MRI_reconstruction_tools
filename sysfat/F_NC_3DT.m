function F = F_NC_3DT(freqs, Ns, Nro, Nt, Nx, Ny, Nz, Nc, varargin)
%function F = F_NC_3DT(freqs, Ns, Nro, Nt, Nx, Ny, Nz, Nc)
% Non-Cartesian version
% freqs can be [Nreadout, Nspokes_per_frame, num_frames] OR list mode [num_samples ordered in time]
% or cell to accomodate different number of spokes per frame! {[Nreadout Nspokes_per_this_frame]}_num_frames
% |
% | 	list_mode
% |		output size: [sum(Ns)*Nz Nc]
% |		freqs size: [sum(Ns) 1] 
% | 		use Ns vector to apportion kspace data to each frame
% |             Ns is per slice!

arg.Ns = Ns;
arg.Nro = Nro;
arg.Nt = Nt;
arg.Nx = Nx;
arg.Ny = Ny;
arg.Nz = Nz;
arg.Nr = arg.Nx * arg.Ny * arg.Nz;
arg.Nc = Nc;
arg.mask = true(arg.Nx, arg.Ny, arg.Nz, arg.Nt, arg.Nc);
arg = vararg_pair(arg, varargin);

for frame_ndx = 1:arg.Nt
	if arg.Ns(frame_ndx) == 0
		A{frame_ndx} = [];
        else
                if frame_ndx == 1
                        curr_ndcs = 1:arg.Ns(1);
                else
                        curr_ndcs = sum(arg.Ns(1:frame_ndx - 1)) + 1 : sum(arg.Ns(1:frame_ndx));
                end
                A{frame_ndx} = GnufftSoS(freqs(curr_ndcs), arg.Ns(frame_ndx), Nx, Ny, Nz);
        end
end
arg.A = A;

% if (arg.Nc > 1)
%         % output is Nt length cell array, each of which is [Ns_f Nc]
%         F = fatrix2('idim', [arg.Nx arg.Ny arg.Nt arg.Nc] ,'arg',arg,'odim', ...
%                 [sum(arg.Ns) arg.Nc], 'forw', @F_NC_3DT_forw, ...
%                 'back', @F_NC_3DT_back, 'imask', arg.mask);
% else
%         F = fatrix2('idim', [arg.Nx arg.Ny arg.Nt] ,'arg',arg,'odim', ...
%                 [sum(arg.Ns) arg.Nc], 'forw', @F_NC_3DT_forw, ...
%                 'back', @F_NC_3DT_back, 'imask', arg.mask);
% end

F = fatrix2('idim', [arg.Nx arg.Ny arg.Nz arg.Nt arg.Nc] ,'arg',arg,'odim', ...
        [sum(arg.Ns)*arg.Nz arg.Nc], 'forw', @F_NC_3DT_forw, ...
        'back', @F_NC_3DT_back, 'imask', arg.mask);


end

% y = G * x
function S = F_NC_3DT_forw(arg, s)

S = [];
for coil_ndx = 1:arg.Nc
        coil_S = [];
        for frame_ndx = 1:arg.Nt
                if ~isempty(arg.A{frame_ndx})
                        curr = col(s(:,:,:, frame_ndx, coil_ndx));
                        curr_S = arg.A{frame_ndx}*(curr);
                        coil_S = [coil_S; col(curr_S)];
                end
        end
        S = [S coil_S];
end

end

% x = G' * y
function s = F_NC_3DT_back(arg, S)
s = zeros(arg.Nx, arg.Ny, arg.Nz, arg.Nt, arg.Nc);
for coil_ndx = 1:arg.Nc
        for frame_ndx = 1:arg.Nt
                if ~isempty(arg.A{frame_ndx})
                        if frame_ndx == 1
                                curr_ndcs = 1:arg.Ns(1);
                        else
                                curr_ndcs = sum(arg.Ns(1:frame_ndx - 1)) + 1 : sum(arg.Ns(1:frame_ndx));
                        end
                        if max(curr_ndcs) > size(S,1)
                                keyboard;
                        end
                        curr_S = S(curr_ndcs, coil_ndx);                        
                        s(:,:,:,frame_ndx, coil_ndx) = reshape(arg.A{frame_ndx}'*curr_S, arg.Nx, arg.Ny, arg.Nz);
                else
                        s(:,:,:,frame_ndx, coil_ndx) = zeros(arg.Nx, arg.Ny, arg.Nz);
                end
	end
end
end
