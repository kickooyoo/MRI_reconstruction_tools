function dcf = radial_dcf(Nro)
%function dcf = radial_dcf(Nro)
%
%

center = 0.5/(Nro/2);
k = linspace(-0.5, 0.5, Nro)';
dcf = abs(k)/max(k) + center;


