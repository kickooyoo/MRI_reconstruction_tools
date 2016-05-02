function F = F_NC_3D(freqs, Ns, Nro, Nt, Nx, Ny, Nc, varargin)
%function F = F_NC_3D(freqs, Ns, Nro, Nt, Nx, Ny, Nc)
% Non-Cartesian version
% freqs can be [Nreadout, Nspokes_per_frame, num_frames] OR list mode [num_samples ordered in time]
% or cell to accomodate different number of spokes per frame! {[Nreadout Nspokes_per_this_frame]}_num_frames
% |
% | Several modes available:
% | 	list_mode
% |		output size: [sum(Ns) Nc]
% |		freqs size: [sum(Ns) 1] 
% | 		use Ns vector to apportion kspace data to each frame
% | 	regular mode 
% | 		output size: [Nro Nspokespf Nt Nc]
% | 		freqs size: [Nro Nspokespf Nt]
% | 		same number of spokes/kspace assigned to each frame
% | 	new mode? how to deal with coronal problem?
% | 
% | 
% | 
% | 
% to do: revert to true 3D + time

arg.Ns = Ns;
arg.Nro = Nro;
arg.Nt = Nt;
arg.Nx = Nx;
arg.Ny = Ny;
arg.Nr = arg.Nx*arg.Ny;
arg.Nc = Nc;
arg.mask = true(arg.Nx, arg.Ny, arg.Nt, arg.Nc);
arg = vararg_pair(arg, varargin);
arg.list_mode = (ndims_ns(freqs) == 1) & ~iscell(freqs);
if ~arg.list_mode %% reconcile! same as Ns
	arg.Nspokespf = size(freqs,2);
end

for frame_ndx = 1:arg.Nt
	if arg.Ns(frame_ndx) == 0
		A{frame_ndx} = 0;
	else
		Jd = [6,6];
		Kd = floor([Nx, Ny]*2);
		if arg.list_mode
			if frame_ndx == 1
				curr_ndcs = 1:arg.Ns(1);
			else
				curr_ndcs = sum(arg.Ns(1:frame_ndx - 1)) + 1 : sum(arg.Ns(1:frame_ndx));
			end
			om = [real(col(freqs(curr_ndcs))), imag(col(freqs(curr_ndcs)))]*2*pi;
		elseif iscell(freqs)
			om = [real(col(freqs{frame_ndx})), imag(col(freqs{frame_ndx}))]*2*pi;
		else
			om = [real(col(freqs(:,:,frame_ndx))), imag(col(freqs(:,:,frame_ndx)))]*2*pi;
		end
		A{frame_ndx} = Gnufft(true(Nx, Ny), {om; [Nx, Ny]; Jd; Kd; [Nx, Ny]/2; 'table'; 2^10; 'minmax:kb'});
	end
end
arg.A = A;

if arg.list_mode
	if (arg.Nc > 1)
		% output is Nt length cell array, each of which is [Ns_f Nc]
		F = fatrix2('idim', [arg.Nx arg.Ny arg.Nt arg.Nc] ,'arg',arg,'odim', ...
			[sum(arg.Ns) arg.Nc], 'forw', @F_NC_3D_forw, ...
			'back', @F_NC_3D_back, 'imask', arg.mask);
	else
		F = fatrix2('idim', [arg.Nx arg.Ny arg.Nt] ,'arg',arg,'odim', ...
			[sum(arg.Ns) arg.Nc], 'forw', @F_NC_3D_forw, ...
			'back', @F_NC_3D_back, 'imask', arg.mask);
	end
elseif iscell(freqs)	
	if (arg.Nc > 1)
		F = fatrix2('idim', [arg.Nx arg.Ny arg.Nt arg.Nc] ,'arg',arg,'odim', ...
			[arg.Nro arg.Nspokespf arg.Nt arg.Nc], 'forw', @F_NC_3D_forw, ...
			'back', @F_NC_3D_back, 'imask', arg.mask);
	else
		F = fatrix2('idim', [arg.Nx arg.Ny arg.Nt] ,'arg',arg,'odim', ...
			[arg.Nro arg.Nspokespf arg.Nt], 'forw', @F_NC_3D_forw, 'back', ...
			@F_NC_3D_back, 'imask', arg.mask);
	end
else
	if (arg.Nc > 1)
		F = fatrix2('idim', [arg.Nx arg.Ny arg.Nt arg.Nc] ,'arg',arg,'odim', ...
			[arg.Nro arg.Nspokespf arg.Nt arg.Nc], 'forw', @F_NC_3D_forw, ...
			'back', @F_NC_3D_back, 'imask', arg.mask);
	else
		F = fatrix2('idim', [arg.Nx arg.Ny arg.Nt] ,'arg',arg,'odim', ...
			[arg.Nro arg.Nspokespf arg.Nt], 'forw', @F_NC_3D_forw, 'back', ...
			@F_NC_3D_back, 'imask', arg.mask);
	end
end

end

% y = G * x
function S = F_NC_3D_forw(arg, s)

if arg.list_mode
	S = [];
else
	S = zeros(arg.Nro,arg.Nspokespf,arg.Nt,arg.Nc);
end
for coil_ndx = 1:arg.Nc
	if arg.list_mode
		coil_S = [];
	end
	for frame_ndx = 1:arg.Nt
		curr = col(s(:,:, frame_ndx, coil_ndx));
		if arg.list_mode
			curr_S = arg.A{frame_ndx}*(curr);
			coil_S = [coil_S; col(curr_S)];
		else
			S(:,:, frame_ndx, coil_ndx) = reshape(arg.A{frame_ndx}*(curr), arg.Nro, arg.Nspokespf);
		end
	end
	if arg.list_mode
		S = [S coil_S];
	end
end

end

% x = G' * y
function s = F_NC_3D_back(arg, S)
% display('YAYYYYY GOT INTO THE TRANSPOSE');
s = zeros(arg.Nx, arg.Ny, arg.Nt, arg.Nc);
for coil_ndx = 1:arg.Nc
	for frame_ndx = 1:arg.Nt
		if arg.list_mode
			if frame_ndx == 1
				curr_ndcs = 1:arg.Ns(1);
			else
				curr_ndcs = sum(arg.Ns(1:frame_ndx - 1)) + 1 : sum(arg.Ns(1:frame_ndx));
			end
			if max(curr_ndcs) > size(S,1)
				keyboard;
			end
			curr_S = S(curr_ndcs, coil_ndx);
		%end
		%if arg.list_mode
			curr = curr_S;%(:, coil_ndx);
		else
			curr = col(S(:,:, frame_ndx, coil_ndx));
		end
		s(:,:, frame_ndx, coil_ndx) = reshape(arg.A{frame_ndx}'*(curr), arg.Nx, arg.Ny);
	end
end
% display('FINISHED TRANSPOSE')
% keyboard;
end
