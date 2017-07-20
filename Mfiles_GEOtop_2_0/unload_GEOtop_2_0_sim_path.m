function unload_GEOtop_2_0_sim_path(mysim)

%function GEOtop_2_0_sim_path(mysim)
% adapted for GEOTOP GEOtop2_0 outputs
% unload the path of the mysimulation mysim
% Copyright, 2009 Giacomo Bertoldi
%
% Inputs:
% 'mysim': string name of the simulation
%
% Outputs:
% clear variable simpath

% load simulation path
simpath=[mysim];
rmpath(genpath(simpath));
clear simpath;


return
