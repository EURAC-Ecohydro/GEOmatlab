function [teta] = teta_from_psi(psi,i,s,r,a,n,m,pmin,Ss)

%[teta] = teta_psi(psi,i,s,r,a,n,m,pmin,Ss)
% Calculates VWC theta [] from suction potential psi [mm] 
% 
% From code GEOtop 1.225-9 From pedo.funct.c
% from function teta_psi
%
% still works only with scalars 
% psi should be negative
%
%   Inputs			
% psi	suction potential [mm]		
% i	volumetric ice content []		
% s	volumetric saturated max water+ice content []		
% r	volumetric residual minimum water content []		
% a	Van Genuchten alpha [mm^-1]		
% n	Van Genuchten n []		
% m	Van Genuchten m=1-1/n		
% pmin	Minimum Psi [mm]		
% Ss	Specific Storativity [mm^-1]		
		
sat=0;

psisat=(((((1.0-i/(s-r))^(-1.0/m))-1.0)^(1.0/n)))*(-1.0/a);

if(psi>psisat) sat=1; end

if(psi<pmin) psi=pmin; end

if(sat==0)
    if(psi>-1.E-6)
        TETA=1.0;
    else
        TETA=1.0/((1.0+((a*(-psi))^n))^m);
    end
    teta=r+TETA*(s-r);
else
    teta= s-i + Ss * (psi-psisat);
end

return 


