In this repository ase saved several Matlab Scripts to import 1D point output from GEOtop 2.0 model and plot it in Matlab.

The main script to load GEOtop timeseries data  is:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 load_GEOtop_2_0.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loads GEOTOP GEOtop 2.0 temporal inputs/outputs
% 
% it parese the geotop.inpts file and then imports the 
% timeseries output files
%
% Copyright, 2017 Giacomo Bertoldi, Stefano della Chiesa, 2011 Ruth
% Mugford
%
% you should edit the variable 'simname' and the 'filroot'
% in this m-file to specify the simulation path
%
% you should also set the appropriate flags

The main script to plot 1D point GEOtop outptut file data is:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_GEOtop_2_0_point.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright, 2010 Giacomo Bertoldi;
%
% plots GEOTOP 2_0 temporal outputs for a chosen control point ;
% load before load_GEOtop_2_0.m
% you need to define the chosen point # in the script

The main script to plot 1D point profiles of soil moisture, soil pressure, soil temperature is:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 plot_GEOtop_2_0_profiles.m
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Copyright, 2010 Giacomo Bertoldi
%
% plot 2_0 GEOTOP point profiles for every layer
% of: soil temparature, sil pressure, soil total water content, soil ice
% content
%
% load before load_GEOtop_2_0


The following functions are called:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function simpath=load_GEOtop_2_0_sim_path(mysim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% adapted x geotop GEOtop2_0 simulations
% load the path of the simulation mysim
% Copyright, 2009 Giacomo Bertoldi
%
% Inputs:
% 'mysim': string name of the simulation
%
% Outputs:
% simpath: simulation path

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function unload_GEOtop_2_0_sim_path(mysim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% adapted for GEOTOP GEOtop2_0 outputs
% unload the path of the mysimulation mysim
% Copyright, 2009 Giacomo Bertoldi
%
% Inputs:
% 'mysim': string name of the simulation
%
% Outputs:
% clear variable simpath

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mydata, mytimes, mydates, mylabels]=import_GEOtop_2_0_data(myfile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% adapted x geotop GEOtop2_0 simulations
% import a GEOtop2_0 time series file
% Copyright, 2009 Giacomo Bertoldi
%
% Inputs:
% 'myfile': string name of the file to import
%
% Outputs:
% mydata: array with numeric data
% mytimes: array with the time in days (Matlab standard date time)
% mydates: string array with the dates
% mylabels: string array with column labels
