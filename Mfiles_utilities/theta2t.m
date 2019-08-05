function t=theta2t(theta,p_atm)
%function t=theta2t(theta,p_atm)
% matlab function to calculate T from T potential
% Copyright, 2009 Giacomo Bertoldi
%
% theta: t potential [C]
% t: t [C]
% p_atm: pressure [hpa]

Tk=273.15;
cp=1005;

t=(theta+Tk)/((1000/p_atm)^(287.04/cp))-Tk;

return 

