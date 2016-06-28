function mag = straighten_complex_sinusoid(signals, varargin)
%function mag = straighten_complex_sinusoid(signals, varargin)
% when signal is looping in complex plane at an angle, abs, real, and imag fail to represent amplitude
% inputs:
%       signals (complex double) [Nt Nsignals]
%
% varargin:
%       tvar_phi (boolean)
%               estimate time-varying rotational angle
%       reg (double) [2 1] 
%               sqrt regularization param, 1st: Ax, Ay, 2nd: phi if
%               tvar_phi
%       figs_on (boolean)
%               default: false
%
% parameterize waveform with:
% s(t) = rot(A(t) * e^(i*omega*t + theta), phi(t))
% [sx(t) sy(t)] = rot([Ax(t)*cos(omega*t + theta), Ay(t)*sin(omega*t +
%               theta), phi(t))
%       omega: sinusoidal frequency
%       theta: initial angle
%       Ax: time varying amplitude in x
%       Ay
%       phi: time varying rotational angle

arg.figs_on = false;
arg.reg = [2 3]*1; % sqrt lambda
arg.tvar_phi = false; % doesn't work great yet
arg = vararg_pair(arg, varargin);

[Nt, Nsignals] = size(signals);

for ii = 1:Nsignals
        
        % center signal, separate x and y
        rsig = real(signals(:,ii));
        rsig = rsig - mean(rsig);
        isig = imag(signals(:,ii));
        isig = isig - mean(isig);
        
        % initialize nonlinear vars
        if arg.tvar_phi
                x = [13; 0; 0.5*ones(Nt, 1); 0.1*ones(Nt, 1); zeros(Nt, 1)];
                % r1^2 *||C Ax||^2 + r1^2 * ||C Ay||^2 + r2^2 * ||C phi||^2
                reg_term = @(x) [arg.reg(1) * (x(4 : 3 + Nt - 1) - circshift(x(4 : 3 + Nt - 1), 1)); %Ax
                        arg.reg(1) * (x(4 + Nt : 3 + 2*Nt - 1) - circshift(x(4 + Nt : 3 + 2*Nt - 1), 1)); %Ay
			arg.reg(1) * (sqrt(abs(x(4 : 3 + Nt - 1)).^2 + abs(x(4 + Nt : 3 + 2*Nt - 1)).^2) ...
			- sqrt(abs(circshift(x(4 : 3 + Nt - 1), 1)).^2 + abs(circshift(x(4 + Nt : 3 + 2*Nt - 1), 1)).^2)); % magnitude
                        arg.reg(2) * (x(4 + 2*Nt : end) - circshift(x(4 + 2*Nt : end), 1))]; %phi
                sinu_diff = @(x) [construct_from_params(x, arg); zeros(4*(Nt - 1),1)] - ...
                        [col(rsig); col(isig); reg_term(x)];
                phi_size = Nt;
        else
                x = [13; 0; 0.5*ones(Nt, 1); 0.1*ones(Nt, 1); 0];
                % r1^2 *||C Ax||^2 + r1^2 * ||C Ay||^2 
                reg_term = @(x) [arg.reg(1) * (x(4 : 3 + Nt - 1) - circshift(x(4 : 3 + Nt - 1), 1)); %Ax
                        arg.reg(1) * (x(4 + Nt : 3 + 2*Nt - 1) - circshift(x(4 + Nt : 3 + 2*Nt - 1), 1)); %Ay;
			arg.reg(1) * (sqrt(abs(x(4 : 3 + Nt - 1)).^2 + abs(x(4 + Nt : 3 + 2*Nt - 1)).^2) ...
			- sqrt(abs(circshift(x(4 : 3 + Nt - 1), 1)).^2 + abs(circshift(x(4 + Nt : 3 + 2*Nt - 1), 1)).^2))]; % magnitude
                sinu_diff = @(x) [construct_from_params(x, arg); zeros(3*(Nt - 1),1)] - ...
                        [col(rsig); col(isig); reg_term(x)];
                phi_size = 1;
        end

        % constrain amplitudes to be positive and greater than margin
        [rpeak_vals, peak_ndcs] = findpeaks(abs(rsig));
        Axmin = min(rpeak_vals) * 0.5;
        Axmax = max(rpeak_vals) * 1.5;
        [ipeak_vals, peak_ndcs] = findpeaks(abs(isig));
        Aymin = min(ipeak_vals) * 0.5;
        Aymax = max(ipeak_vals) * 1.5;
        lb = [0; 0; Axmin*ones(Nt,1); Aymin*ones(Nt, 1); zeros(phi_size,1)];
        ub = [Inf; 2*pi; Axmax*ones(Nt, 1); Aymax*ones(Nt, 1); 2*pi*ones(phi_size, 1)];
        
        nlvar = lsqnonlin(sinu_diff, x, lb, ub);
        
        t = (0:Nt - 1).';
        omega_hat = nlvar(1);
        theta_hat = nlvar(2);
        Ax_hat = nlvar(3 : 3 + Nt - 1);
        Ay_hat = nlvar(3 + Nt : 3 + 2*Nt - 1);
        if arg.tvar_phi
                phi_hat = nlvar(3 + 2*Nt : end);
        else
                phi_hat = nlvar(end);
        end
        
        x_mag = Ax_hat .* cos(omega_hat*t + theta_hat);
        i_mag = Ay_hat .* sin(omega_hat*t + theta_hat);
	
	rot = zeros(2*Nt);
	if ~arg.tvar_phi
		phi_hat = repmat(phi_hat, [Nt 1]);
	end
	for tt = 1:Nt
		rot(2*(tt - 1) + 1: tt*2, 2*(tt - 1) + 1 : 2*tt) = [cos(phi_hat(tt)) -sin(phi_hat(tt)); sin(phi_hat(tt)) cos(phi_hat(tt))];
	end
	slanted_mag = rot*[x_mag; i_mag];
        
	if sum(Ax_hat > Ay_hat) > Nt/2
                mag(:,ii) = x_mag;
        else
                mag(:,ii) = i_mag;
        end
        if arg.figs_on
                figure; subplot(1,2,1); plot(mag(:,ii));
                subplot(1,2,2); plot(signals(:,ii));
                hold on; plot(x_mag + 1i*i_mag, 'r');
		hold on; plot(slanted_mag(1:Nt) + 1i*slanted_mag(Nt + 1:end), 'g');
                axis equal             
        end
end
keyboard
end


function slanted_sinu = construct_from_params(x, arg)
%omega: sinusoidal frequency
%theta: initial angle
%Ax: time varying amplitude in x
%Ay
%phi: time varying rotational angle

if arg.tvar_phi
        Nt = (length(x) - 2)/3;
else
        Nt = (length(x) - 3)/2;
end
if mod(Nt, 1) ~= 0, keyboard; end
omega = x(1);
theta = x(2);
Ax = x(3 : 3 + Nt - 1);
Ay = x(3 + Nt : 3 + 2*Nt - 1);
if arg.tvar_phi
        phi = x(3 + 2*Nt : end);
else
        phi = repmat(x(end), [Nt 1]);
end

% rotation matrix
rot = zeros(2*Nt);
for tt = 1:Nt
	rot(2*(tt - 1) + 1: tt*2, 2*(tt - 1) + 1 : 2*tt) = [cos(phi(tt)) -sin(phi(tt)); sin(phi(tt)) cos(phi(tt))];
end

t = (0:Nt - 1).';
stretch_sinu = [Ax .* cos(omega*t + theta); Ay .* sin(omega*t + theta)];
slanted_sinu = rot*stretch_sinu;
end
