function rh=q2rh(qa,p_atm,Ta)
%function rh=q2rh(q,p_atm)
% matlab function to calculate rh from q
% Copyright, 2009 Giacomo Bertoldi
%
% rh: relative humidity [%]
% qa: q [g/kg]
% p_atm: pressure [hpa]
% Ta: air tempearture [C]

Tk=273.15;

a1=373.15;
a2=13.3185;
a3=3.952;
a4=1.9335;
a5=0.5196;

b1=1013.25;
b2=13.3185;
b3=1.9760;
b4=0.6445;
b5=0.1299;

Tak=Ta+Tk;      
tr=1-373.15./Tak;
estar=b1.*exp(b2.*tr - b3.*tr.^2 - b4.*tr.^3 - b5.*tr.^4);
qasat=0.622.*estar./(p_atm-0.378.*estar);
rh=qa./1000./qasat.*100;



return 

