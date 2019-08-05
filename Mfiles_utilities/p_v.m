
function ee=p_v(Qa,p_atm)
%function pv=p_v(Qa,p_atm)
% matlab function to calculate air vapour pressure [hpa]
% Copyright, 2010 Giacomo Bertoldi
%
% ee: air vapour pressure [hpa]
% Qa: air specific humidity [kg/kg]
% p_atm: pressure [hpa]


ee=Qa*p_atm/(0.622+Qa*0.378);

return 



% GEOtop 1.1

%ee=Qa*Ppoint/(0.622+Qa*0.378);