function F = GnufftSoS(freqs, Ns, Nx, Ny, Nz, varargin)
% function F = GnufftSoS(freqs, Ns, Nx, Ny, Nz, varargin)
%
% freqs: row vec
% real: kx, imag: ky
%
% one frame, one coil
% 2D freqs only!
% since stack of stars is same in z, use same 2D Gnufft

arg.Ns = Ns;
arg.Nx = Nx;
arg.Ny = Ny;
arg.Nz = Nz;
arg.mask = true(arg.Nx, arg.Ny, arg.Nz);
arg.threeD = true;
arg = vararg_pair(arg, varargin);
freqs = col(freqs);

if arg.threeD
	assert(mod(Nz, 2) == 0, '3D option in GnufftSoS only available for even Nslices (Nz)');
	Nfreq = numel(freqs);
	kz = col(double(repmat(permute(col(-Nz/2:Nz/2-1)/Nz, [2 1]), [Nfreq 1])));
	k3D = repmat(freqs, [1 Nz]);%exp(1i*kz).*repmat(freqs, [1 Nz]);
	%om_z = double(repmat(permute(col(-Nz/2:Nz/2-1)/Nz, [2 3 1]), [Nro Nspokes 1]));
	kx = real(col(k3D));
	ky = imag(col(k3D));
	om = 2*pi*[kx, ky, kz];
	simple_dims = [Nx Ny Nz];
	Jd = [6 6 6];
else
	om = [real(col(freqs)), imag(col(freqs))]*2*pi;
	simple_dims = [Nx, Ny];
	Jd = [6 6];
end
arg.A = Gnufft(true(simple_dims), {om; simple_dims; Jd; ceil(simple_dims*1.5); simple_dims/2; 'table'; 2^10; 'minmax:kb'});

F = fatrix2('idim', [arg.Nx arg.Ny arg.Nz] ,'arg',arg,'odim', ...
        [arg.Ns*arg.Nz 1], 'forw', @GnufftSoS_forw, ...
        'back', @GnufftSoS_back, 'imask', arg.mask);

end

function S = GnufftSoS_forw(arg, s)

if arg.threeD
	S = arg.A*s;
else
	S = [];
	for z = 1:arg.Nz
		curr = col(s(:,:,z));
		curr_S = arg.A*curr;
		S = [S; col(curr_S)];
	end
end

end

function s = GnufftSoS_back(arg, S)
if arg.threeD
	s = arg.A'*S;
else
	for z = 1:arg.Nz
		curr_ndcs = (z-1)*arg.Ns + (1:arg.Ns);
		s(:,:,z) = arg.A'*S(curr_ndcs);
	end
end
end

