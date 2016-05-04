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
% | 
% | if stack-of-stars (default), freqs and Ns should not include Nz

if (nargin == 1) && strcmp(freqs, 'test')
	F = F_NC_3DT_test();
	return;
end

arg.Ns = Ns; % does not include Nz!
arg.Nro = Nro;
arg.Nt = Nt;
arg.Nx = Nx;
arg.Ny = Ny;
arg.Nz = Nz;
arg.Nr = arg.Nx * arg.Ny * arg.Nz;
arg.Nc = Nc;
arg.mask = true(arg.Nx, arg.Ny, arg.Nz, arg.Nt, arg.Nc);
arg.SoS = true; % stack of stars
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
		if arg.SoS
			A{frame_ndx} = GnufftSoS(freqs(curr_ndcs), arg.Ns(frame_ndx), Nx, Ny, Nz);
		else
			keyboard
			% A{frame_ndx} = nufft?? % to do
		end
        end
end
arg.A = A;

if (arg.Nc > 1)
	F = fatrix2('idim', [arg.Nx arg.Ny arg.Nz arg.Nt arg.Nc] ,'arg',arg,'odim', ...
		[sum(arg.Ns)*arg.Nz arg.Nc], 'forw', @F_NC_3DT_forw, ...
		'back', @F_NC_3DT_back, 'imask', arg.mask);
else
	if (arg.Nt > 1)
		F = fatrix2('idim', [arg.Nx arg.Ny arg.Nz arg.Nt] ,'arg',arg,'odim', ...
			[sum(arg.Ns)*arg.Nz 1], 'forw', @F_NC_3DT_forw, ...
			'back', @F_NC_3DT_back, 'imask', arg.mask);
	else
		F = fatrix2('idim', [arg.Nx arg.Ny arg.Nz] ,'arg',arg,'odim', ...
			[sum(arg.Ns)*arg.Nz 1], 'forw', @F_NC_3DT_forw, ...
			'back', @F_NC_3DT_back, 'imask', arg.mask);
	end
end

%F = fatrix2('idim', [arg.Nx arg.Ny arg.Nz arg.Nt arg.Nc] ,'arg',arg,'odim', ...
%        [sum(arg.Ns)*arg.Nz arg.Nc], 'forw', @F_NC_3DT_forw, ...
%        'back', @F_NC_3DT_back, 'imask', arg.mask);


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
                                curr_ndcs = 1:arg.Ns(1)*arg.Nz;
                        else
                                curr_ndcs = sum(arg.Ns(1:frame_ndx - 1)*arg.Nz) + 1 : sum(arg.Ns(1:frame_ndx)*arg.Nz);
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

function F = F_NC_3DT_test()

Nx = 64;
Ny = 50;
Nz = 30;
img = phantom3d('Modified Shepp-Logan', max([Nx Ny Nz]));
img = crop_outside(img, [Nx Ny Nz]); 
[xx, yy, zz] = ndgrid(1:Nx, 1:Ny, 1:Nz);
ph = 0.2*xx+0.5*yy-0.1*zz;
ph = ph/max(col(ph))*0.9*pi;
img = img.*exp(1i*ph);

Nro = round(Nx*1);
Nspokes = 200;
Nslice = Nz;
Nslice_PF = round(0.8*Nslice);
Nslice_nocollect = Nslice - Nslice_PF;
grad_shift = 0; % 0.75?
k = get_GA_coords(Nro, Nspokes, 1, 'grad_shift', grad_shift);

Ns = Nro*Nspokes;
Nt = 1;
Nc = 1;
F = F_NC_3DT(k, Ns, Nro, Nt, Nx, Ny, Nz, Nc);
test_adjoint(F, 'big', 1, 'complex', 1, 'nrep', 10)

return
rawdat = F*img;

rawdat = reshape(rawdat, Nro, Nspokes, Nz);
rawdat_PF = rawdat(:, :, end - Nslice_PF + 1: end);

%% use 3D PF
[~, reconk_3D] = PF_3D(rawdat_PF, [Nro Nspokes Nz], 'PF_location', [0 0 1], 'window_step3', 3);
reconimg_3D = F'*reconk_3D(:);


end
