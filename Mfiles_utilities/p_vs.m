
function ees=p_vs(Ta,p_atm)
%function ees=p_vs(Ta,p_atm)
% matlab function to calculate air vapour pressure at saturation [hPa]
% Copyright, 2010 Giacomo Bertoldi
%
% ees: air vapour pressure at saturation [hpa]
% Ta: air temperature [C]
% p_atm: pressure [hpa]

A=6.1121*(1.0007+3.46E-6*p_atm);
b=17.502;
c=240.97;
ees=A.*exp(b.*Ta./(c+Ta));

return 

% GEOtop 1.1
% 
% void sat_vap_pressure(double *p, double *dp_dT, double T, double P){	//water vapour pressure p(T,P)
% //pressure in [mbar] - p vapour pressure - P atmospheric pressure, temperature in [C]
% 	double A, b, c;
% 	A=6.1121*(1.0007+3.46E-6*P);
% 	b=17.502;
% 	c=240.97;
% 	*p=A*exp(b*T/(c+T));
% 	*dp_dT=*p*(b/(c+T)-b*T/pow(c+T,2.0));
% }