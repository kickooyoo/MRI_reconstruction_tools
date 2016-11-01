function FS = F_NC_S_5D(freqs, sense_maps, Ns, Nro, Nt, Nresp, varargin)
% function FS = F_NC_S_5D(freqs, sense_maps, Ns, Nro, Nt, Nresp, varargin)
% Non-Cartesian version
% freqs can be [Nreadout, Nspokes_per_frame, num_frames] OR list mode [num_samples ordered in time]
% or cell to accomodate different number of spokes per frame! {[Nreadout Nspokes_per_this_frame]}_num_frames
% | inputs:
% | 	freqs ([1 Nresp] cell of [Nro*Nspokes(resp)*Nslice] complex double)
% |     sense_maps [Nx Ny Nz Nc]
% |	Ns ([1 Nresp] cell of [1 Nt] int)
% |		DOES THIS INCLUDE NRO AND NSLICE FACTOR??
% | varargin:
% | 	list_mode
% |		output size: [sum(Ns) Nc]
% |		freqs size: Nresp length cell, each [sum(Ns_resp) 1],
% |                     Ns_resp per slice!
% | 		use Ns vector to apportion kspace data to each frame
% |	sampling (bool) [Nro Nspokes Nslice Nresp]
% |		if empty, use pre-binned data as described above
% | 		if provided (xsparse)
% | 			expect freqs [Nro Nspokes]
% | 			leave Ns empty
% |	doboth (bool) [1 2]
% |		do F, do S, default: true(1,2)
% | 	small_mask (bool) [Nx Ny]
% | output of fatrix is single column of samples
% | 	ordered:
% | 		[readout, spoke, slice, coil]
% | 
arg.Ns = Ns;
arg.Nro = Nro;
arg.Nt = Nt;
arg.Nx = size(sense_maps, 1); % are these set if doboth(2) == 0?
arg.Ny = size(sense_maps, 2);
arg.Nz = size(sense_maps, 3);
arg.Nc = size(sense_maps, 4);
arg.Nresp = Nresp;
arg.Nr = arg.Nx * arg.Ny * arg.Nz;
arg.sampling = [];
arg.smaps = sense_maps;
arg.doboth = [true true];
arg.small_mask = []; % save memory
arg.verbose = false;
arg.pf = true;
arg.Nworkers = 8;
arg = vararg_pair(arg, varargin);

if(arg.pf)
	pool = gcp('nocreate');
	if numel(pool) == 0
		pool = parpool();
	end 
	if isempty(arg.Nworkers)
		arg.Nworkers = pool.NumWorkers;
	end
end

if ~isempty(arg.sampling) % apply sampling for user
	[freqs, arg.Ns] = apply_sampling(freqs, arg); 
end

all_Ns = zeros(Nt, Nresp);


if arg.pf
	[A, arg.all_Ns] = construct_all_nufft_pf(freqs, arg);
	arg.cum_Ns = cumsum(col(arg.all_Ns));
else
	[A, arg.all_Ns] = construct_all_nufft(freqs, arg);
	arg.cum_Ns = reshape(cumsum(col(arg.all_Ns)), arg.Nt, arg.Nresp);
end
arg.A = A;

if arg.doboth(1)
	odims = [sum(col(arg.all_Ns))*arg.Nz arg.Nc];
else
	odims = [arg.Nx arg.Ny arg.Nz arg.Nt arg.Nresp arg.Nc];
end
if arg.doboth(2)
	idims = [arg.Nx arg.Ny arg.Nz arg.Nt arg.Nresp];
else
	idims = [arg.Nx arg.Ny arg.Nz arg.Nt arg.Nresp arg.Nc];
end
if ~isempty(arg.small_mask) && ((ndims(arg.small_mask) ~= 2) || ( ~all(size(arg.small_mask) == idims(1:2))))
	display('dims of mask not right');
	keyboard;
end
if ~isempty(arg.small_mask)
	idims = [sum(col(arg.small_mask)) idims(3:end)];
end

if arg.pf
	FS = fatrix2('idim', idims, 'arg', arg,'odim', odims, 'forw', @F_NC_S_5D_forw_pf, ...
		'back', @F_NC_S_5D_back_pf);
else
	FS = fatrix2('idim', idims, 'arg', arg,'odim', odims, 'forw', @F_NC_S_5D_forw, ...
		'back', @F_NC_S_5D_back);
end

end

% y = FSx
function y = F_NC_S_5D_forw(arg, x)

if ~isempty(arg.small_mask)
	x = embed(x, arg.small_mask);
end
y = [];
for coil_ndx = 1:arg.Nc
	if arg.verbose, tic, end
        coil_S = [];
        for resp_ndx = 1:arg.Nresp
		for frame_ndx = 1:arg.Nt
			if ~isempty(arg.A{frame_ndx, resp_ndx})
				if arg.doboth(2)
					curr_s = x(:,:,:, frame_ndx, resp_ndx) .* arg.smaps(:,:,:, coil_ndx);
				else
					curr_s = x(:,:,:, frame_ndx, resp_ndx, coil_ndx);
				end
				if arg.doboth(1)
					curr_S = arg.A{frame_ndx, resp_ndx}*col(curr_s);
				else
					curr_S = col(curr_s);
				end
				coil_S = [coil_S; col(curr_S)];
			end
		end
	end
	if arg.verbose, display(sprintf('done with %d/%d coils in %d sec', coil_ndx, arg.Nc, toc)), end
        y = [y coil_S];
end
end

function y = F_NC_S_5D_forw_pf(arg, x)

if ~isempty(arg.small_mask)
	x = embed(x, arg.small_mask);
end

% parfor will not slice fields of structs!
A = arg.A;
if arg.doboth(2)
	x = reshape(x, arg.Nx, arg.Ny, arg.Nz, arg.Nt*arg.Nresp);
else
	x = reshape(x, arg.Nx, arg.Ny, arg.Nz, arg.Nt*arg.Nresp, arg.Nc);
end

% don't want to duplicate smpas to each worker...
% Nc Nx Ny Nz Nw 
% but if I explicitly construct Sx first, then I have ...
% Nx Ny Nz Nt Nr Nc
% Nw vs. Nt*Nr

% concern: need coil as 2nd dim of y
parfor (frame_resp_ndx = 1:arg.Nresp*arg.Nt, arg.Nworkers)
	if ~isempty(A{frame_resp_ndx})
		if arg.doboth(2)
			curr_s = repmat(squeeze(x(:,:,:,frame_resp_ndx)), [1 1 1 arg.Nc]).*arg.smaps;
		else
			curr_s = x(:,:,:, frame_resp_ndx, :);
		end
		if arg.doboth(1)
			% reshape to use does_many over Nz, Nc
			curr_s = A{frame_resp_ndx}*reshape(curr_s, arg.Nx, arg.Ny, arg.Nz*arg.Nc);
			% curr_s  [arg.Ns{resp_ndx}(frame_ndx) Nz*Nc]
			% where arg.Nx{resp_ndx}(frame_ndx) is Nro*Nspokes of frame, resp
			curr_s = reshape(curr_s, size(curr_s, 1)*arg.Nz, arg.Nc);
		end
		curr_S_pf{frame_resp_ndx} = curr_s; 			
	end
end

% construct y because had to save slices over frame_resp_ndx
y = [];
for frame_resp_ndx = 1:arg.Nresp*arg.Nt
	y = [y; curr_S_pf{frame_resp_ndx}];
end

end


% x = F'S'y
function x = F_NC_S_5D_back(arg, y)

if arg.doboth(2)
	x = zeros(arg.Nx, arg.Ny, arg.Nz, arg.Nt, arg.Nresp, 'single');
else
	x = zeros(arg.Nx, arg.Ny, arg.Nz, arg.Nt, arg.Nresp, arg.Nc, 'single');
end

for resp_ndx = 1:arg.Nresp
        for frame_ndx = 1:arg.Nt
		if arg.verbose, tic, end
		if ~isempty(arg.A{frame_ndx, resp_ndx})
			if arg.doboth(1)
				if (frame_ndx == 1) && (resp_ndx == 1)
					curr_ndcs = 1:arg.Nz*arg.cum_Ns(frame_ndx, resp_ndx);
				elseif frame_ndx == 1
					curr_ndcs = arg.Nz*arg.cum_Ns(arg.Nt, resp_ndx - 1) + 1: arg.Nz*arg.cum_Ns(frame_ndx, resp_ndx);
				else
					curr_ndcs = arg.Nz*arg.cum_Ns(frame_ndx - 1, resp_ndx) + 1: arg.Nz*arg.cum_Ns(frame_ndx, resp_ndx);
				end
				if max(curr_ndcs) > size(y,1)
					keyboard;
				end
				curr_S = y(curr_ndcs, :); 
				curr_A = arg.A{frame_ndx, resp_ndx};
				if ~all(size(curr_S) == curr_A.odim)
					keyboard
				end
				% induce 'does_many' over Nz, Nc
				curr_s = arg.A{frame_resp_ndx}'*reshape(curr_S, size(curr_S,1)/arg.Nz, arg.Nz*arg.Nc);
				% output [Nx Ny Nz*Nc]
			else
				curr_s =  y(:,:,:, frame_ndx, resp_ndx,:);
			end
			small_s = reshape(curr_s, arg.Nx, arg.Ny, arg.Nz, arg.Nc);
		else
			small_s = zeros(arg.Nx, arg.Ny, arg.Nz, arg.Nc);
		end
		if arg.doboth(2)
			small_prod = conj(arg.smaps) .* small_s;
			x(:,:,:, frame_ndx, resp_ndx) = sum(small_prod, 4);
		else
			x(:,:,:, frame_ndx, resp_ndx, :) = permute(small_s, [1 2 3 5 6 4]);
		end
	end
	if arg.verbose, display(sprintf('done with resp %d/%d in %d sec', resp_ndx, arg.Nresp, toc)), end
end
if ~isempty(arg.small_mask)
	x = masker(x, arg.small_mask);
end

end


function x = F_NC_S_5D_back_pf(arg, y)

% parfor will not slice fields of structs!
A = arg.A;
if ~arg.doboth(1)
	y = reshape(y, arg.Nx, arg.Ny, arg.Nz, arg.Nt*arg.Nresp, arg.Nc);
end
keyboard
if 1
parfor (frame_resp_ndx = 1:arg.Nresp*arg.Nt, 2)%arg.Nworkers)
	if ~isempty(A{frame_resp_ndx})
		if arg.doboth(1)
			if (frame_resp_ndx == 1)
				curr_ndcs = 1:arg.Nz*arg.cum_Ns(frame_resp_ndx);
			else
				curr_ndcs = arg.Nz*arg.cum_Ns(frame_resp_ndx - 1) + 1: arg.Nz*arg.cum_Ns(frame_resp_ndx);
			end
			if max(curr_ndcs) > size(y,1)
				keyboard;
			end
			curr_S = y(curr_ndcs, :);
			% induce 'does_many' over Nz, Nc
			curr_s = A{frame_resp_ndx}'*reshape(curr_S, size(curr_S,1)/arg.Nz, arg.Nz*arg.Nc);
			% output [Nx Ny Nz*Nc]
		else
			curr_s = y(:,:,:, frame_resp_ndx,:);
		end
		small_s = reshape(curr_s, arg.Nx, arg.Ny, arg.Nz, arg.Nc);
	else
		small_s = zeros(arg.Nx, arg.Ny, arg.Nz, arg.Nc);
	end
	small_s_pf(:,:,:,frame_resp_ndx,:) = small_s;
end
else
	small_s_pf = zeros(arg.Nx, arg.Ny, arg.Nz, arg.Nt*arg.Nresp, arg.Nc);
end
if arg.doboth(2)
	% for loop to avoid duplicating memory hog smaps
	for frame_resp_ndx = 1:arg.Nresp*arg.Nt
		x(:,:,:,frame_resp_ndx) = sum(conj(arg.smaps) .* squeeze(small_s_pf(:,:,:,frame_resp_ndx,:)),4);
	end
	x = reshape(x, arg.Nx, arg.Ny, arg.Nz, arg.Nt, arg.Nresp);
else	
	x = reshape(small_s_pf, arg.Nx, arg.Ny, arg.Nz, arg.Nt, arg.Nresp, arg.Nc);
end
if ~isempty(arg.small_mask)
	x = masker(x, arg.small_mask);
end

end

function [A, all_Ns] = construct_all_nufft(freqs, arg)

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
				%A{frame_ndx, resp_ndx} = GnufftSoS(curr_freqs(curr_ndcs), curr_Ns(frame_ndx), arg.Nx, arg.Ny, arg.Nz);
				% default: 'table', 2^10, 'minmax:kb'
				dims = [arg.Nx arg.Ny];
				k = curr_freqs(curr_ndcs);
				om = [real(col(k)) imag(col(k))]*2*pi;
				Jd = [6 6];
				% rely on does_many for Nz and Nc
				A{frame_resp, resp_ndx} = Gnufft({om, dims, Jd, ceil(dims*1.5), max(floor(dims/2),1), 'table', 2^10, 'minmax:kb'});
			end
		end
	end
	arg.all_Ns = all_Ns;

end


function [A, all_Ns] = construct_all_nufft_pf(freqs, arg)

	for frame_resp_ndx = 1:arg.Nresp*arg.Nt
		resp_ndx = ceil(frame_resp_ndx / arg.Nt);
		frame_ndx = mod(frame_resp_ndx - 1, arg.Nt) + 1;

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
		if curr_Ns(frame_ndx) == 0
			A{frame_resp_ndx} = [];
		else
			A{frame_resp_ndx} = 1;
			if frame_ndx == 1
				curr_ndcs = 1:curr_Ns(1);
			else
				curr_ndcs = sum(curr_Ns(1:frame_ndx - 1)) + 1 : sum(curr_Ns(1:frame_ndx));
			end
			%A{frame_resp_ndx} = GnufftSoS(curr_freqs(curr_ndcs), curr_Ns(frame_ndx), arg.Nx, arg.Ny, arg.Nz);
			% default: 'table', 2^10, 'minmax:kb'
			dims = [arg.Nx arg.Ny];
			k = curr_freqs(curr_ndcs);
			om = [real(col(k)) imag(col(k))]*2*pi;
			Jd = [6 6];
                        % rely on 'does_many' option for Nz and Nc
			A{frame_resp_ndx} = Gnufft({om, dims, Jd, ceil(dims*1.5), max(floor(dims/2),1), 'table', 2^10, 'minmax:kb'});
		end
	end
	arg.all_Ns = all_Ns;

end

function [freqs, Ns] = apply_sampling(freqs, arg)

% do I need to know Nspokes for this?
% move sampling out of arg? hijack an input var?

if issparse(arg.sampling)
	if size(arg.sampling,1) ~= size(arg.sampling,2)
		display('non-square sampling sparse matrix, why?');
		keyboard
	end
	Nspokes = size(arg.sampling,1) / (arg.Nro * arg.Nz * arg.Nresp);
	M = Nspokes/arg.Nt;
	if mod(M, 1) ~= 0
		display('to do: non-integer factor between Nt and Nspokes');
		keyboard;
	end
	spoke_ndcs = (1:Nspokes);
	Nsparse_dresp = arg.Nro * Nspokes * arg.Nz;
	for ii = 1:arg.Nresp
		resp_ndcs = Nsparse_dresp * (ii - 1) + (1:Nsparse_dresp);
		% [Nro Nslice] number of frames associated with given respiratory state
		in_resp = reshape(full(diag(arg.sampling(resp_ndcs, resp_ndcs))), arg.Nro, arg.Nz, Nspokes); 
		Nsamp_per_spoke = Nro*squeeze(in_resp(1,1,:));
		Ns{ii} = sum(reshape(Nsamp_per_spoke, M, arg.Nt), 1); 
		freqs_per_resp{ii} = col(freqs(:, spoke_ndcs(logical(in_resp(1,1,:)))));
	end
	freqs = freqs_per_resp;
else	
	Nspokes = size(arg.sampling, 2);
	if ~all([arg.Nro, Nspokes, arg.Nz, arg.Nresp] == [size(arg.sampling, 1) size(arg.sampling, 2) size(arg.sampling, 3), size(arg.sampling, 4)])
		display('size mismastch');
		keyboard
	end
	M = Nspokes/arg.Nt;
	if mod(M, 1) ~= 0
		display('to do: non-integer factor between Nt and Nspokes');
		keyboard;
	end
	spoke_ndcs = (1:Nspokes);
	for ii = 1:arg.Nresp
		in_resp = squeeze(sum(arg.sampling(:,:,:,ii),2)); % [Nro Nslice] number of frames associated with given respiratory state
		Nsamp_per_spoke = Nro*squeeze(arg.sampling(1,:,1,ii));
		Ns{ii} = sum(reshape(Nsamp_per_spoke, M, arg.Nt), 1); 
		freqs_per_resp{ii} = col(freqs(:, spoke_ndcs(logical(arg.sampling(1,:,1,ii)))));
	end
	freqs = freqs_per_resp;
end
end
