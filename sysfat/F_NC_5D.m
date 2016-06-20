function F = F_NC_5D(freqs, Ns, Nro, Nt, Nresp, Nx, Ny, Nz, Nc, varargin)
%function F = F_NC_5D(freqs, Ns, Nro, Nt, Nresp, Nx, Ny, Nz, Nc, varargin)
% Non-Cartesian version
% freqs can be [Nreadout, Nspokes_per_frame, num_frames] OR list mode [num_samples ordered in time]
% or cell to accomodate different number of spokes per frame! {[Nreadout Nspokes_per_this_frame]}_num_frames
% | inputs:
% | 	freqs ([1 Nresp] cell of [Nro*Nspokes(resp)*Nslice] complex double)
% |	Ns ([1 Nresp] cell of [1 Nt] int)
% |		DOES THIS INCLUDE NRO AND NSLICE FACTOR??
% |
% | 	list_mode
% |		output size: [sum(Ns) Nc]
% |		freqs size: Nresp length cell, each [sum(Ns_resp) 1],
% |                     Ns_resp per slice!
% | 		use Ns vector to apportion kspace data to each frame
% |	sampling (bool) [Nro Nspokes Nslice Nresp]
% |		if empty, use pre-binned data as described above
% | 		if provided (sparse)
% | 			expect freqs [Nro Nspokes]
% | 			leave Ns empty
% | 
arg.Ns = Ns;
arg.Nro = Nro;
arg.Nt = Nt;
arg.Nx = Nx;
arg.Ny = Ny;
arg.Nz = Nz;
arg.Nresp = Nresp;
arg.Nr = arg.Nx * arg.Ny * arg.Nz;
arg.Nc = Nc;
arg.sampling = [];
% arg.mask = true(arg.Nx, arg.Ny, arg.Nz, arg.Nt, arg.Nr, arg.Nc);
arg = vararg_pair(arg, varargin);

% do I need to know Nspokes for this?
% move sampling out of arg? hijack an input var?
if ~isempty(arg.sampling) % apply sampling for user
	if issparse(arg.sampling)
		if size(arg.sampling,1) ~= size(arg.sampling,2)
			display('non-square sampling sparse matrix, why?');
			keyboard
		end
		Nspokes = size(arg.sampling,1) / (Nro * Nz * Nresp);
		M = Nspokes/Nt;
		if mod(M, 1) ~= 0
			display('to do: non-integer factor between Nt and Nspokes');
			keyboard;
		end
		spoke_ndcs = (1:Nspokes);
		Nsparse_dresp = Nro * Nspokes * Nz;
		for ii = 1:Nresp
			resp_ndcs = Nsparse_dresp * (ii - 1) + (1:Nsparse_dresp);
			in_resp = reshape(full(diag(arg.sampling(resp_ndcs, resp_ndcs))), Nro, Nz, Nspokes); % [Nro Nslice] number of frames associated with given respiratory state
			Nsamp_per_spoke = Nro*squeeze(in_resp(1,1,:));
			arg.Ns{ii} = sum(reshape(Nsamp_per_spoke, M, Nt), 1); 
			freqs_per_resp{ii} = col(freqs(:, spoke_ndcs(logical(in_resp(1,1,:)))));
		end
		freqs = freqs_per_resp;
	else	
		Nspokes = size(arg.sampling, 2);
		if ~all([Nro, Nspokes, Nz, Nresp] == [size(arg.sampling, 1) size(arg.sampling, 2) size(arg.sampling, 3), size(arg.sampling, 4)])
				display('size mismastch');
			keyboard
		end
		M = Nspokes/Nt;
		if mod(M, 1) ~= 0
			display('to do: non-integer factor between Nt and Nspokes');
			keyboard;
		end
		spoke_ndcs = (1:Nspokes);
		for ii = 1:Nresp
			in_resp = squeeze(sum(arg.sampling(:,:,:,ii),2)); % [Nro Nslice] number of frames associated with given respiratory state
			Nsamp_per_spoke = Nro*squeeze(arg.sampling(1,:,1,ii));
			arg.Ns{ii} = sum(reshape(Nsamp_per_spoke, M, Nt), 1); 
			freqs_per_resp{ii} = col(freqs(:, spoke_ndcs(logical(arg.sampling(1,:,1,ii)))));
		end
		freqs = freqs_per_resp;
	end
end

all_Ns = zeros(Nt, Nresp);
for resp_ndx = 1:arg.Nresp
        curr_Ns = arg.Ns{resp_ndx};
	if length(curr_Ns) ~= arg.Nt
		display('bad input Ns');
		keyboard
	end
        all_Ns(:,resp_ndx) = curr_Ns;
        curr_freqs = freqs{resp_ndx};
	if numel(curr_freqs) ~= sum(curr_Ns)
		display('k and Ns mismatch');
		keyboard
	end
        for frame_ndx = 1:arg.Nt
                if curr_Ns(frame_ndx) == 0
                        A{frame_ndx, resp_ndx} = [];
                else
                        A{frame_ndx, resp_ndx} = 1;
                        if frame_ndx == 1
                                curr_ndcs = 1:curr_Ns(1);
                        else
                                curr_ndcs = sum(curr_Ns(1:frame_ndx - 1)) + 1 : sum(curr_Ns(1:frame_ndx));
                        end
			% default: 'table', 2^10, 'minmax:kb'
                        A{frame_ndx, resp_ndx} = GnufftSoS(curr_freqs(curr_ndcs), curr_Ns(frame_ndx), Nx, Ny, Nz);
                end
        end
end
arg.A = A;
arg.all_Ns = all_Ns;
arg.cum_Ns = reshape(cumsum(col(arg.all_Ns)), arg.Nt, arg.Nresp);

F = fatrix2('idim', [arg.Nx arg.Ny arg.Nz arg.Nt arg.Nresp arg.Nc] ,'arg',arg,'odim', ...
        [sum(col(arg.all_Ns))*arg.Nz arg.Nc], 'forw', @F_NC_5D_forw, ...
        'back', @F_NC_5D_back);%, 'imask', arg.mask);


end

% y = G * x
function S = F_NC_5D_forw(arg, s)

S = [];
for coil_ndx = 1:arg.Nc
        coil_S = [];
        for resp_ndx = 1:arg.Nresp
                for frame_ndx = 1:arg.Nt
                        if ~isempty(arg.A{frame_ndx, resp_ndx})
                                curr_s = col(s(:,:,:, frame_ndx, resp_ndx, coil_ndx));
                                curr_S = arg.A{frame_ndx, resp_ndx}*curr_s;
                                coil_S = [coil_S; col(curr_S)];
                        end
                end
        end
        S = [S coil_S];
end

end

% x = G' * y
function s = F_NC_5D_back(arg, S)
s = zeros(arg.Nx, arg.Ny, arg.Nz, arg.Nt, arg.Nresp, arg.Nc);
for coil_ndx = 1:arg.Nc
        for resp_ndx = 1:arg.Nresp
                for frame_ndx = 1:arg.Nt
                        if ~isempty(arg.A{frame_ndx, resp_ndx})
                                if (frame_ndx == 1) && (resp_ndx == 1)
                                        curr_ndcs = 1:arg.Nz*arg.cum_Ns(frame_ndx, resp_ndx);
                                elseif frame_ndx == 1
                                        curr_ndcs = arg.Nz*arg.cum_Ns(arg.Nt, resp_ndx - 1) + 1: arg.Nz*arg.cum_Ns(frame_ndx, resp_ndx);
                                else
                                        try
                                        curr_ndcs = arg.Nz*arg.cum_Ns(frame_ndx - 1, resp_ndx) + 1: arg.Nz*arg.cum_Ns(frame_ndx, resp_ndx);
                                        catch
                                                keyboard
                                        end
                                end
                                if max(curr_ndcs) > size(S,1)
                                        keyboard;
                                end
                                curr_S = S(curr_ndcs, coil_ndx);
                                curr_A = arg.A{frame_ndx, resp_ndx};
                                if ~all(size(curr_S) == curr_A.odim)
                                        keyboard
                                end
                                curr_s = arg.A{frame_ndx, resp_ndx}'*curr_S;
                                s(:,:,:, frame_ndx, resp_ndx, coil_ndx) = reshape(curr_s, arg.Nx, arg.Ny, arg.Nz);
                        else
                                s(:,:,:, frame_ndx, resp_ndx, coil_ndx) = zeros(arg.Nx, arg.Ny, arg.Nz);
                        end
                end
        end
end
end
