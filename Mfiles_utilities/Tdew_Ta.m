function Td=Tdew_Ta(Ta,Rh,p_atm)
%function Td=Tdew_Ta(Ta,Rh,p_atm)
% matlab function to calculate dew point temperature
% To calculate Tdew if you know Rh[], p_atm[hPa], Ta[C]
% Copyright, 2010 Giacomo Bertoldi
%
% Td: dew point temperature [C]
% Ta: air temperature [C]
% Rh: air realtive humidity (fraction!) []
% p_atm: pressure [hpa]
%

ees=p_vs(Ta,p_atm); 
Qas=Q_as(ees,p_atm); 
Qa=Q_a(Rh,Qas);
pv=p_v(Qa,p_atm);
Td=Tdew(pv,p_atm);

return