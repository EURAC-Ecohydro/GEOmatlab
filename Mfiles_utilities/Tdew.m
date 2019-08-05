function Td=Tdew(p_v,p_atm)
%function Td=Tdew(p_v,p_atm)
% matlab function to calculate dew point temperature
% Copyright, 2010 Giacomo Bertoldi
%
% Td: dew point temperature [C]
% p_v: air vapour pressure [hpa]
% p_atm: pressure [hpa]
%
% To calculate Tdew if you know Rh[], p_atm[hPa], Ta[C]
% ees=p_vs(Ta,p_atm), 
% Qas=Q_as(ees,p_atm), 
% Qa=Q_a(Rh,Qas), 
% pv=p_v(Qa,p_atm), 
% dew=Tdew(pv,p_atm)


A=6.1121*(1.0007+3.46E-6*p_atm);
b=17.502;
c=240.97;
Td=c.*log(p_v./A)/(b-log(p_v./A));

return 

% GEOtop 1.1

%void sat_vap_pressure_inv(double *T, double p, double P){	//temperature(vapour pressure p,P)
% 	double A, b, c;
% 	A=6.1121*(1.0007+3.46E-6*P);
% 	b=17.502;
% 	c=240.97;
% 	*T=c*log(p/A)/(b-log(p/A));
% }

% //calculate dew temperature (otherwise replace Tdew with met->Tgrid) to distinguish between rain and snow
% 				sat_vap_pressure(&ee,&dee,Tpoint,Ppoint);	
% 				Qa=RHpoint*spec_humidity(ee, Ppoint);			
% 				ee=Qa*Ppoint/(0.622+Qa*0.378);
% 				sat_vap_pressure_inv(&Tdew,ee,Ppoint);

