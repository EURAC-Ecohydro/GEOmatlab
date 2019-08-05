function P=Pz(Po,Zo,Z)
%P=Pz(Po,Zo,Z)
% matlab function to calculate standard pressure elevation dependance
% Copyright, 2009 Giacomo Bertoldi
%
% Po: pressure at Z=Zo (hPa)
% Zo: input elevation m
% Z: output elevation m
% P: pressure at Z=Z (hPa)

P=Po*exp(-(Z-Zo)*0.00013);



return 

