function ET_mm=Et_conversion(LE,Ta,dt)
%function ET_mm=Et_conversion(LE,Ta,dt)
% matlab function to get ET in mm from LE in W/m2
% Copyright, 2011 Giacomo Bertoldi
%
% LE: latent heat flux [W/m2]
% Ta: air tempearture [C]
% dt: time step in second over which is averaged LEs

% Latent heat of evapotranspiration L [J/Kg]
L=2501000+(2406000-2501000)*(Ta)/40;

rho_w=1000; % water density kg/m3

ET_mm = 1000 * LE./(L.*rho_w)*dt;


return 

