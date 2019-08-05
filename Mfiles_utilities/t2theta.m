function theta=t2theta(t,p_atm)
%function theta=t2theta(t,p_atm)
% matlab function to calculate T potential from T 
% Copyright, 2009 Giacomo Bertoldi
%
% theta: t potential [C]
% t: t [C]
% p_atm: pressure [hpa]

Tk=273.15;
cp=1005;

theta=(t+Tk)*((1000/p_atm)^(287.04/cp))-Tk;

return 

