function simpath=load_GEOtop_2_0_sim_path(mysim)

%function simpath=load_GEOtop2_0_path(mysim)
% adapted x geotop GEOtop2_1 simulations
% load the path of the simulation mysim
% Copyright, 2009 Giacomo Bertoldi
%
% Inputs:
% 'mysim': string name of the simulation
%
% Outputs:
% simpath: simulation path

% get the current path
mypath=pwd;

% load simulation path
simpath=mysim;
addpath(genpath(simpath));



return
