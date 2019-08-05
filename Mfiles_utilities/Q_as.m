
function Qas=Q_as(ees,p_atm)
%function Qas=Q_as(ees,p_atm)
% matlab function to calculate air specific humidity at saturation [kg/kg]
% Copyright, 2010 Giacomo Bertoldi
%
% ees: air vapour pressure at saturation [hpa]
% Qas: air specific humidity at saturation [kg/kg]
% p_atm: pressure [hpa]


Qas=0.622*ees/(p_atm-0.378*ees);

return 

% GEOtop 1.1

% double spec_humidity(double p, double P){
% 	double Q;
% 	Q=0.622*p/(P-0.378*p);
% 	return(Q);
% }