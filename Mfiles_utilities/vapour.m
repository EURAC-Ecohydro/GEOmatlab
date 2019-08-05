function [Pv dPvdT]=vapour(Ta,Pa)
%function [Pv dPvdT]=vapour(Ta,Pa)
% matlab function to calculate the water vapour saturation pressure
% Copyright, 2009 Giacomo Bertoldi
%
% Ta: T air  [C]
% Pa: atmospheric pressure [hpa]
% Pv: water vapour saturation pressure [hpa]
% dPv_dT: dPv/dTa [hpa/K]

A=6.1121*(1.0007+3.46E-6*Pa);
b=17.502;
c=240.97;

Pv=A*exp(b*Ta/(c+Ta));
dPvdT=Pv*(b/(c+Ta)-b*Ta/(c+Ta)^2);

return
