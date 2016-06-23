function mag = straighten_complex_sinusoidal(signal)
%function mag = straighten_complex_sinusoidal(signal)
% when signal is looping in complex plane at an angle, abs, real, and imag fail to represent amplitude

Nt = length(signal);
rsig = real(signal);
isig = imag(signal);
x = [13; 0; 0.5*ones(Nt, 1); 0.1*ones(Nt, 1); zeros(Nt, 1)];
%test = construct_from_params(x);

sinu_diff = @(x) construct_from_params(x) - [col(rsig); col(isig)];%signal;
nlvar = lsqnonlin(sinu_diff, x);

t = (0:Nt - 1).';
omega_hat = nlvar(1);
theta_hat = nlvar(2);
Ax_hat = nlvar(3 : 3 + Nt - 1);
Ay_hat = nlvar(3 + Nt : 3 + 2*Nt - 1);
phi_hat = nlvar(3 + 2*Nt : end);

if sum(Ax_hat > Ay_hat) > Nt/2
	mag = Ax_hat .* cos(omega_hat*t + phi_hat);
else
	mag = Ay_hat .* sin(omega_hat*t + phi_hat);
end

keyboard;

end


function slanted_sinu = construct_from_params(x)
%omega, theta, Ax, Ay, phi)

Nt = (length(x) - 2)/3;
if mod(Nt, 1) ~= 0, keyboard; end
omega = x(1);
theta = x(2);
Ax = x(3 : 3 + Nt - 1);
Ay = x(3 + Nt : 3 + 2*Nt - 1);
phi = x(3 + 2*Nt : end);

rot = zeros(2*Nt);
for tt = 1:Nt
	rot(2*(tt - 1) + 1: tt*2, 2*(tt - 1) + 1 : 2*tt) = [cos(phi(tt)) -sin(phi(tt)); sin(phi(tt)) cos(phi(tt))];
end
t = (0:Nt - 1).';
stretch_sinu = [Ax .* cos(omega*t + phi); Ay .* sin(omega*t + phi)];
stack_sinu = rot*stretch_sinu;

slanted_sinu = stack_sinu;%stack_sinu(1:Nt) + 1i*(stack_sinu(Nt + 1:end));
end
