function Qa=Q_a(Rh,Qas)
%function Qa=Q_a(Rh,Qas)
% matlab function to calculate air specific humidity [kg/kg]
% Copyright, 2010 Giacomo Bertoldi
%
% Rh: air relative humidity []
% Qas: air specific humidity at saturation [kg/kg]


Qa=Rh.*Qas;

return 


% GEOtop 1.1

%Qa=RHpoint*spec_humidity(ee, Ppoint);	