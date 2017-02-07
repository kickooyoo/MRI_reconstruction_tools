function FS = F_NC_S_5D(freqs, sense_maps, Ns, Nro, Nt, Nresp, varargin)
% function FS = F_NC_S_5D(freqs, sense_maps, Ns, Nro, Nt, Nresp, varargin)
% Non-Cartesian in kx-ky, Cartesian in kz
% 
% freqs can be [Nreadout, Nspokes_per_frame, num_frames] OR list mode [num_samples ordered in time]
% or cell to accomodate different number of spokes per frame! {[Nreadout Nspokes_per_this_frame]}_num_frames
%
% | inputs:
% | 	freqs [Nro Nspokes(sorted byresp)]  assuming same across Cartesian Nz, (complex double)
% |     sense_maps [Nx Ny Nz Nc]
% |	Ns [Nt Nresp] (int)
% |		number of slices associated with Nt, Nresp (not including Nro, Nslice factor)
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
% | input of fatrix is 
% |	 	[Nx Ny Nz Nt Nresp]
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
arg.small_mask = []; % for vectorized, masked x
arg.verbose = false;
arg.pf = true;
arg.Nworkers = [];
arg.small_imask = []; % implicit masker on image x
arg.debug = 0;
arg.chat = 1;
arg = vararg_pair(arg, varargin);

if(arg.pf)
	pool = gcp('nocreate');
	if numel(pool) == 0
		if isempty(arg.Nworkers)
			pool = parpool();
			arg.Nworkers = pool.NumWorkers;
		else
			pool = parpool(arg.Nworkers);
		end
	end 
end

if ~isempty(arg.sampling) % apply sampling for user
	[freqs, arg.Ns] = apply_sampling(freqs, arg); 
end

if arg.pf
	[A, arg.cum_Ns] = construct_all_nufft_pf(freqs, arg);
else
	[A, arg.cum_Ns] = construct_all_nufft(freqs, arg);
end
arg.A = A;

if arg.doboth(1)
	odims = [arg.Nro arg.cum_Ns(end) arg.Nz arg.Nc];
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

% ---------------------------------------------------------------------------------------------------------
% y = FSx
% x: [Nx Ny Nz Nt Nresp]
% Sx: [Nx Ny Nz Nt Nresp Nc]
% y: [Nro Nspokes Nslice Nc]
function y = F_NC_S_5D_forw(arg, x) %----------------------------------------
if ~isempty(arg.small_mask)
	x = embed(x, arg.small_mask);
end
%if ~isempty(arg.small_imask)
%	x = masker(x, arg.small_imask);
%	smaps = masker(arg.smaps, arg.small_imask);
%else
%	x = reshape(x, arg.Nx*arg.Ny, arg.Nz, arg.Nt, arg.Nresp);
%	smaps = reshape(arg.smaps, arg.Nx*arg.Ny, arg.Nz, arg.Nc); 
%end
if arg.doboth(1)
	y_dims = [arg.Nro arg.cum_Ns(end), arg.Nz, arg.Nc];
else
	y_dims = [arg.Nx arg.Ny arg.Nz arg.Nt arg.Nresp arg.Nc];
end
y = zeros(y_dims, 'single');
for coil_ndx = 1:arg.Nc
	if arg.verbose, tic, end
        for resp_ndx = 1:arg.Nresp
		for frame_ndx = 1:arg.Nt
			if ~isempty(arg.A{frame_ndx, resp_ndx})
				if arg.doboth(2)
					curr_s = x(:,:,:, frame_ndx, resp_ndx) .* arg.smaps(:,:,:, coil_ndx);
				else
					curr_s = x(:,:,:, frame_ndx, resp_ndx, coil_ndx);
				end
				if arg.doboth(1)
					% need to duce does_many NOT WORKING
					curr_S = arg.A{frame_ndx, resp_ndx}*curr_s;
					if arg.debug && all([coil_ndx frame_ndx resp_ndx] == [1 1 1])
						tmp = reshape(arg.A{frame_ndx, resp_ndx}'*curr_S, size(curr_s));
						subplot(1,3,1); im('mid3', curr_s);
						subplot(1,3,2); im('mid3', tmp);
						drawnow;
						%keyboard
					end
					spoke_ndcs = get_spoke_ndcs([frame_ndx resp_ndx], arg.cum_Ns, arg.Nt);
					y(:,spoke_ndcs,:,coil_ndx) = reshape(curr_S, [arg.Nro numel(spoke_ndcs) arg.Nz]);
				else
					y(:,:,:, frame_ndx, resp_ndx, coil_ndx) = curr_s;
				end
			end
		end
	end
	if arg.verbose, display(sprintf('done with %d/%d coils in %d sec', coil_ndx, arg.Nc, toc)), end
end
end

% ---------------------------------------------------------------------------------------------------------
function y = F_NC_S_5D_forw_pf(arg, x) % ------------------------------------------

if arg.chat, display('beginning FS forw pf'), end
if ~isempty(arg.small_mask)
	x = embed(x, arg.small_mask);
end

% parfor will not slice fields of structs!
A = arg.A;
Nx = arg.Nx;
Ny = arg.Ny;
Nz = arg.Nz;
Nc = arg.Nc;
Nt = arg.Nt;
Nresp = arg.Nresp;
doboth = arg.doboth;
smaps = arg.smaps;

if arg.doboth(2)
	x = reshape(x, arg.Nx, arg.Ny, arg.Nz, arg.Nt*arg.Nresp);
else
	x = reshape(x, arg.Nx, arg.Ny, arg.Nz, arg.Nt*arg.Nresp, arg.Nc);
end

% don't want to duplicate smaps to each worker...
% Nc Nx Ny Nz Nw 
% but if I explicitly construct Sx first, then I have ...
% Nx Ny Nz Nt Nr Nc
% Nw vs. Nt*Nr
% outside for loop over coil

% concern: need coil as 2nd dim of y

if doboth(2)
	parfor (coil_ndx = 1:arg.Nc, arg.Nworkers)
		if arg.chat, display(sprintf('entered %d/%d thread in FS forw pf, S', coil_ndx, arg.Nc)), end
		if mod(coil_ndx, 4) == 0, whos, end
		s(:,:,:,:,coil_ndx) = x.*repmat(smaps(:,:,:,coil_ndx), [1 1 1 Nt*Nresp]);
	end
else
	s = x;
end

parfor (frame_resp_ndx = 1:arg.Nresp*arg.Nt, arg.Nworkers)
	if arg.chat, display(sprintf('entered %d/%d thread in FS forw pf', frame_resp_ndx, arg.Nresp*arg.Nt)), end
	if ~isempty(A{frame_resp_ndx})
		if doboth(1)
			% reshape to use does_many over Nz, Nc
			curr_S = A{frame_resp_ndx}*reshape(s(:,:,:, frame_resp_ndx, :), Nx, Ny, Nz*Nc);
			if arg.chat, display(sprintf('done applying F in %d/%d thread in FS forw pf', frame_resp_ndx, arg.Nresp*arg.Nt)), end
			% curr_S  [arg.Ns{resp_ndx}(frame_ndx) Nz*Nc]
			% where arg.Nx{resp_ndx}(frame_ndx) is Nro*Nspokes of frame, resp
			curr_S = reshape(curr_S, size(curr_S, 1)*Nz, Nc);
		else
			curr_S = s(:,:,:, frame_resp_ndx, :);
		end
		curr_S_pf{frame_resp_ndx} = curr_S; 			
	else
		curr_S_pf{frame_resp_ndx} = [];
	end
end

% construct y because had to save slices over frame_resp_ndx
if arg.doboth(1)
	for frame_resp_ndx = 1:arg.Nresp*arg.Nt
			spoke_ndcs = get_spoke_ndcs(frame_resp_ndx, arg.cum_Ns, arg.Nt);
			y_dims = [arg.Nro numel(spoke_ndcs), arg.Nz, arg.Nc];
			y(:,spoke_ndcs,:,:) = reshape(curr_S_pf{frame_resp_ndx}, y_dims);
	end
else
	y_dims = [arg.Nx arg.Ny arg.Nz arg.Nc];
	for frame_resp_ndx = 1:arg.Nresp*arg.Nt
		y(:,:,:,frame_resp_ndx,:) =  reshape(curr_S_pf{frame_resp_ndx}, y_dims);
	end
	y = reshape(y, arg.Nx, arg.Ny, arg.Nz, arg.Nt, arg.Nresp, arg.Nc);
end

end


% ---------------------------------------------------------------------------------------------------------
% x = F'S'y
% y: [Nro Nspokes Nslice Nc]
% S'y: [Nro Nspokes Nslice]
% x: [Nx Ny Nz Nt Nresp]
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
				spoke_ndcs = get_spoke_ndcs([frame_ndx resp_ndx], arg.cum_Ns, arg.Nt);
				if max(spoke_ndcs) > arg.cum_Ns(end) % == Nspokes
					keyboard;
				end
				curr_S = y(:, spoke_ndcs, :, :); 
				% induce 'does_many' over Nz, Nc
				curr_s = arg.A{frame_ndx, resp_ndx}'*reshape(curr_S, arg.Nro*arg.Ns(frame_ndx, resp_ndx), arg.Nz*arg.Nc);
				% output [Nx Ny Nz*Nc]
			else
				curr_s =  y(:,:,:, frame_ndx, resp_ndx,:);
			end
			small_s = reshape(curr_s, arg.Nx, arg.Ny, arg.Nz, arg.Nc);
			if arg.debug && all([frame_ndx resp_ndx] == [1 1])
				subplot(1,3,3); im('mid3', small_s(:,:,:,1));
				drawnow;
				keyboard
			end
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


% ---------------------------------------------------------------------------------------------------------
function x = F_NC_S_5D_back_pf(arg, y)

% parfor will not slice fields of structs!
A = arg.A;
cum_Ns = arg.cum_Ns;
Nro = arg.Nro;
Nt = arg.Nt;
Nx = arg.Nx;
Ny = arg.Ny;
Nz = arg.Nz;
Nc = arg.Nc;
doboth = arg.doboth;
if ~arg.doboth(1)
	y = reshape(y, arg.Nx, arg.Ny, arg.Nz, arg.Nt*arg.Nresp, arg.Nc);
end
parfor (frame_resp_ndx = 1:arg.Nresp*arg.Nt, arg.Nworkers)
	if ~isempty(A{frame_resp_ndx})
		if doboth(1)
			spoke_ndcs = get_spoke_ndcs(frame_resp_ndx, cum_Ns, Nt);
			curr_S = y(:, spoke_ndcs, :, :);
			% induce 'does_many' over Nz, Nc
			curr_s = A{frame_resp_ndx}'*reshape(curr_S, Nro*numel(spoke_ndcs), Nz*Nc);
			% output [Nx Ny Nz*Nc]
		else
			curr_s = y(:,:,:, frame_resp_ndx,:);
		end
		small_s = reshape(curr_s, arg.Nx, arg.Ny, arg.Nz, arg.Nc);
	else
		small_s = zeros(Nx, Ny, Nz, Nc);
	end
	small_s_pf(:,:,:,frame_resp_ndx,:) = small_s;
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

% ---------------------------------------------------------------------------------------------------------
function ndcs = get_spoke_ndcs(frame_resp_ndx, cum_Ns, Nt)
if length(frame_resp_ndx) == 1 % frame_resp_ndx
	if (frame_resp_ndx == 1)
		ndcs = 1:cum_Ns(frame_resp_ndx);
	else
		ndcs = cum_Ns(frame_resp_ndx - 1) + 1: cum_Ns(frame_resp_ndx);
	end
elseif length(frame_resp_ndx) == 2 % frame_ndx, resp_ndx
	frame_ndx = frame_resp_ndx(1);
	resp_ndx = frame_resp_ndx(2);
	if (frame_ndx == 1) && (resp_ndx == 1)
		ndcs = 1:cum_Ns(frame_ndx, resp_ndx);
	elseif frame_ndx == 1
		ndcs = cum_Ns(Nt, resp_ndx - 1) + 1: cum_Ns(frame_ndx, resp_ndx);
	else
		ndcs = cum_Ns(frame_ndx - 1, resp_ndx) + 1: cum_Ns(frame_ndx, resp_ndx);
	end
else
	display('invalid use of get_spoke_ndcs');
	keyboard
end

end

% ---------------------------------------------------------------------------------------------------------
function [A, cum_Ns] = construct_all_nufft(freqs, arg)
	cum_Ns = reshape(cumsum(col(arg.Ns)), arg.Nt, arg.Nresp);
	for resp_ndx = 1:arg.Nresp
		for frame_ndx = 1:arg.Nt
			if arg.Ns(frame_ndx, resp_ndx) == 0
				A{frame_ndx, resp_ndx} = [];
			else
				A{frame_ndx, resp_ndx} = 1;
				curr_spokes = get_spoke_ndcs([frame_ndx resp_ndx], cum_Ns, arg.Nt);
				% default: 'table', 2^10, 'minmax:kb'
				dims = [arg.Nx arg.Ny];
				k = freqs(:,curr_spokes);
				om = [real(col(k)) imag(col(k))]*2*pi;
				Jd = [6 6];
				% rely on does_many for Nz and Nc
				A{frame_ndx, resp_ndx} = Gnufft(true(arg.Nx, arg.Ny), {om, dims, Jd, ceil(dims*1.5), max(floor(dims/2),1), 'table', 2^10, 'minmax:kb'});
			end
		end
	end
end

% ---------------------------------------------------------------------------------------------------------

function [A, cum_Ns] = construct_all_nufft_pf(freqs, arg)
	cum_Ns = cumsum(col(arg.Ns));
	for frame_resp_ndx = 1:arg.Nresp*arg.Nt
		if arg.Ns(frame_resp_ndx) == 0
			A{frame_resp_ndx} = [];
		else
			A{frame_resp_ndx} = 1;
			curr_spokes = get_spoke_ndcs(frame_resp_ndx, cum_Ns, arg.Nt);
			% default: 'table', 2^10, 'minmax:kb'
			dims = [arg.Nx arg.Ny];
			k = freqs(:, curr_spokes);
			om = [real(col(k)) imag(col(k))]*2*pi;
			Jd = [6 6];
                        % rely on 'does_many' option for Nz and Nc
			if ~isempty(arg.small_imask)
				A{frame_resp_ndx} = Gnufft(true(arg.Nx, arg.Ny), {om, dims, Jd, ceil(dims*1.5), max(floor(dims/2),1), 'table', 2^10, 'minmax:kb', 'imask', arg.small_imask});
			else
				A{frame_resp_ndx} = Gnufft(true(arg.Nx, arg.Ny), {om, dims, Jd, ceil(dims*1.5), max(floor(dims/2),1), 'table', 2^10, 'minmax:kb'});
			end
		end
	end
end

% ---------------------------------------------------------------------------------------------------------
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
		Ns(:,ii) = sum(reshape(Nsamp_per_spoke, M, arg.Nt), 1); 
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
		Ns(:,ii) = sum(reshape(Nsamp_per_spoke, M, arg.Nt), 1); 
		freqs_per_resp{ii} = col(freqs(:, spoke_ndcs(logical(arg.sampling(1,:,1,ii)))));
	end
	freqs = freqs_per_resp;
end
end
