function k = create_GA_radial_samples(Nspokes, Nro, varargin)
% function k = create_GA_radial_samples(Nspokes, Nro, grad_shift)
% 
% Inputs:
%	Nspokes:	number of spokes
%	Nro:		number of samples in each readout/spoke
%
% Varargin:
%	grad_shift:	shift from gradient delays, default 0
%				[delta between readout samples, 1/Nro]
%	GA_start:	first angle in GA series [radians]
%	figs_on		plots k-space locations
%
% Output:
%	k:			[Nro Nspokes]
%				complex-valued k-space locations for Nspokes*Nro samples 
%				(normalized freq)
%
arg.grad_shift = 0;
arg.GA_start = pi/2;
arg.figs_on = false;
arg = vararg_pair(arg, varargin);

% Calculate angles for Golden-Angle mode
GA = pi*(sqrt(5)-1)/2; % radians
phi = [arg.GA_start:GA:GA*Nspokes];
phi = mod(phi,2*pi);
assert(mod(Nro,2) == 0, 'Number along read out not even!')
delta_ro = 1/Nro;
rho = col([-(Nro/2 - 1):1:Nro/2])*delta_ro;

%apply gradient delay shift
rho = rho + arg.grad_shift*delta_ro;

k = double(rho*exp(-1i*phi));

if arg.figs_on
	figure; plot(k);
	axis tight;
end

