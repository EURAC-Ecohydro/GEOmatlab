% load_GEOtop_2_0.m
%
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
% yuou should also set appropriate paths in the script: GEOMatlab_mfiles_path
%
% you should also set the appropriate flags

%% add the GEOtop_1_2 mfiles path
%clear all
%close all
%clc
%cd('/Users/giacomo/Documents/Matlab/Mfiles_GEOtop_1_2/');
GEOMatlab_mfiles_path
%%
% simulation name

% Win
simname='Matsch_P2_Ref_001'
fileroot='C:\Users\GBertoldi\Documents\Simulations_local\Montacini_elisa\1D\'


% Mac
%simname='Mazia_dstr_WG1_020bis_light'
%fileroot='/Users/giacomo/Documents/Simulations/mazia_dstr_1_2/'


mysim=[fileroot,simname,'/'];

%%

% you need to define the GEOtop version

%version=1.2 % GEOtop1_2 Endrizzi April 2011 version
version=1.225 % GEOtop1_2225 Endrizzi March 2012 version
version=2.0 % GEOtop_2.0 GIT https://github.com/geotopmodel/geotop/tree/geotop20_veg 20.7.2017

if(version==1.225)
    % snowlayer_max number of snow layers in variable Ns
    % in 1.225 version max number fixed
    Ns=11
elseif(version==2.0)
    Ns=10
end


% flag to activate if you have a point simulation
ifpoint=1
% flag to activate if you have miltiple points output files
ifpoints=0
% flag to activate if you have a basin file
ifbasin=1
% flag to activate if you have a Tz file
ifTz=1
% flag to activate if you have a psiz file
ifpsiz=1
% flag to activate if you have a psiztot file
ifpsiztot=0
% flag to activate if you have a flows file
ifflows=1
% flag to activate if you have a thetaz file
ifthetaz=1
% flag to activate if you have a thetaice  file
ifthetaicez=1
% flag to activate if you have a snow file
ifsnow=0
% flag to activate if you would like to calculate daily means
daily=0
% flag to calculate saturation
ifsatz=1

%% Parameters to change above this line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of column to offset for date
%noffset=2; % 1_32 version
% noffset=3; %1.1. version
%noffset=5; %1.145 version
noffset=5; %1.2 1.225 2.0 versions

% load simulation path
load_GEOtop_2_0_sim_path(mysim);

%% define input file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GEOtop 2_0 file names
flow_name=[mysim 'tabs/discharge.txt'];
basin_name=[mysim 'tabs/basin.txt'];
inpts_name=[mysim 'geotop.inpts'];

soil_root=[mysim 'soil/soil00'];

point_root=[mysim 'tabs/point00'];
Tz_root=[mysim 'tabs/soilTz00'];
psiz_root=[mysim 'tabs/psiz00'];
thetaz_root=[mysim 'tabs/thetaliq00'];
thetaicez_root=[mysim 'tabs/thetaice00'];
%snowz_root=[mysim 'tabs/snow00'];


%% open inputs files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp 'base parameters from geotop.inpts'

% read the geotop.inpts file and saves the data in the variable wordc
fid=fopen(inpts_name);
f=fscanf(fid,'%c');
[m n]=size(f);

newline=sprintf('\n');
tab=sprintf('\t');

inc=0;
wrdinc=1;
for chr=1:n
    inc=inc+1;
    if (((f(inc)==' ')) || (f(inc)==newline) || (f(inc)==tab));
        if (inc>=2)
            if (f(inc-1)==',')
                wordc{wrdinc}=f(1:inc-2);
            else
                wordc{wrdinc}=f(1:inc-1);
            end
        else
            wordc{wrdinc}=f(1:inc-1);
        end
        f(1:inc)=[];
        wrdinc=wrdinc+1;
        inc=0;
    end
    
end
[m n]=size(wordc);
fclose(fid);


% read parameters from GEOtop.inpts

% BUG: no more than one tab charachter in GEOtop.inpts !!

%1)Dt'
searchterm='TimeStepEnergyAndWater' ;
for i=1:n
    
    if strcmp(wordc{i},searchterm)==1;
        %         Dt= wordc{i+2};
        Dt= str2double(wordc{i+2});
    end
end

% '2)InitDateDDMMYYYYhhmm' in InitDatehhmm
searchterm='InitDateDDMMYYYYhhmm' ;
for i=1:n
    if strcmp(wordc{i},searchterm)==1;
        InitDateDDMMYYYYhhmm= wordc{i+2};
        InitDatehhmm= wordc{i+3};
    end
end

% '3)EndDateDDMMYYYYhhmm' in EndDatehhmm
searchterm='EndDateDDMMYYYYhhmm' ;
for i=1:n
    if strcmp(wordc{i},searchterm)==1;
        EndDateDDMMYYYYhhmm= wordc{i+2};
        EndDatehhmm= wordc{i+3};
    end
end

%4) dt average for flows ouptut (h) in Dt_out_disc
searchterm='DtPlotDischarge' ;
for i=1:n
    m1=strtrim(wordc{i});
    if strcmp(wordc{i},searchterm)==1;
        Dt_out_disc= str2double(wordc{i+2});
    end
end

% 5) dt average for  pixel ouptut (h) in Dt_out_pix
searchterm='DtPlotPoint' ;
for i=1:n
    if strcmp(wordc{i},searchterm)==1;
        Dt_out_pix= str2double(wordc{i+2});
    end
end

%  dt average for basin outptut (h) in Dt_out_bas
searchterm='DtPlotBasin' ;
for i=1:n
    if strcmp(wordc{i},searchterm)==1;
        Dt_out_bas= str2double(wordc{i+2});
    end
end

%  dt average for OutputSnowMaps (h) in OutputSnowMaps
searchterm='OutputSnowMaps' ;
for i=1:n
    if strcmp(wordc{i},searchterm)==1;
        Dt_out_maps= str2double(wordc{i+2});
    end
end

% FlagSkyViewFactor in FlagSkyViewFactor
searchterm='FlagSkyViewFactor' ;
for i=1:n
    if strcmp(wordc{i},searchterm)==1;
        FlagSkyViewFactor= str2double(wordc{i+2});
    end
end

% PointSim 1 if point simulation, 0 otherwise in ifpoint
searchterm='PointSim' ;
for i=1:n
    if strcmp(wordc{i},searchterm)==1;
        ifpoint= str2double(wordc{i+2})
    end
end

% Some land cover parameters
searchterm='LSAI' ;
for i=1:n
    if strcmp(wordc{i},searchterm)==1;
        LSAI= str2double(wordc{i+2})
    end
end
searchterm='CanopyFraction' ;
for i=1:n
    if strcmp(wordc{i},searchterm)==1;
        CanopyFraction= str2double(wordc{i+2})
    end
end

searchterm='VegHeight' ;
for i=1:n
    if strcmp(wordc{i},searchterm)==1;
        VegHeight= str2double(wordc{i+2})
    end
end

% number of soil types
searchterm='SoilLayerTypes' ;
for i=1:n
    if strcmp(wordc{i},searchterm)==1;
        n_soils= str2double(wordc{i+2})
    end
end

% if there is a soil fle
searchterm='HeaderSoilDz' ;
ifsoil=0;
for i=1:n
    if strcmp(wordc{i},searchterm)==1;
        ifsoil=1;
    end
end



% load number of point outputs files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bug: to read an array you need to separate elements in .inpts
% only with commas, not spaces in elements
searchterm='CoordinatePointX' ;
for i=1:n
    if strcmp(wordc{i},searchterm)==1;
        tmp1= wordc{i+2};
        tmp2=strrep(tmp1, ',', ';');
        CoordinatePointX=str2num(tmp2);
    end
end
searchterm='CoordinatePointY' ;
for i=1:n
    if strcmp(wordc{i},searchterm)==1;
        tmp1= wordc{i+2};
        tmp2=strrep(tmp1, ',', ';');
        CoordinatePointY=str2num(tmp2);
        clear 'tmp1' 'tmp2'
    end
end
npoints=size(CoordinatePointX,1);


searchterm='MeteoStationCoordinateX' ;
for i=1:n
    if strcmp(wordc{i},searchterm)==1;
        tmp1= wordc{i+2};
        tmp2=strrep(tmp1, ',', ';');
        MeteoStationCoordinateX=str2num(tmp2);
    end
end
searchterm='MeteoStationCoordinateY' ;
for i=1:n
    if strcmp(wordc{i},searchterm)==1;
        tmp1= wordc{i+2};
        tmp2=strrep(tmp1, ',', ';');
        MeteoStationCoordinateY=str2num(tmp2);
        clear 'tmp1' 'tmp2'
    end
end
searchterm='MeteoStationElevation' ;
for i=1:n
    if strcmp(wordc{i},searchterm)==1;
        tmp1= wordc{i+2};
        tmp2=strrep(tmp1, ',', ';');
        MeteoStationElevation=str2num(tmp2);
        clear 'tmp1' 'tmp2'
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load other parameters from geotop.inpts

disp 'load snow parameters'

% snowlayer_max number of snow layers in variable Ns

% in 1.225 version max number fixed
if(version~=1.225)
    if(version~=2.0)
        searchterm='NumMaxSnowLayers' ;
        for i=1:n
            m1=strtrim(wordc{i});
            if strcmp(wordc{i},searchterm)==1;
                Ns= str2double(wordc{i+2});
            end
        end
    end
end

% number of meteo stations nstations
searchterm='NumberOfMeteoStations' ;
for i=1:n
    m1=strtrim(wordc{i});
    if strcmp(wordc{i},searchterm)==1;
        nstations= str2double(wordc{i+2});
    end
end

% BUG: those keywords are no allways specifice in geotop.inpts
% thex can be read from DEM or specified in  separate listpoint fine
% the following three readings are not general

% elevation of point output
searchterm='PointElevation' ;
for i=1:n
    m1=strtrim(wordc{i});
    if strcmp(wordc{i},searchterm)==1;
        PointElevation= str2double(wordc{i+2});
    end
end

% slope of point output
PointSlope=-9999;
searchterm='PointSlope' ;
for i=1:n
    m1=strtrim(wordc{i});
    if strcmp(wordc{i},searchterm)==1;
        PointSlope= str2double(wordc{i+2});
    end
end

% aspect of point output
PointAspect=-9999;
searchterm='PointAspect' ;
for i=1:n
    m1=strtrim(wordc{i});
    if strcmp(wordc{i},searchterm)==1;
        PointAspect= str2double(wordc{i+2});
    end
end


%% open flows file in flows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (ifflows)
    disp 'load discharge output in flows'
    
    
    [flows, flows_t, flows_date, flows_labels]=import_GEOtop_2_0_data(flow_name);
    flows_t=flows(:,2);
    
    % quick load outputs
    %[Q, Q_t, Q_date, Q_labels]=import_GEO11_data('/Users/giacomo/Documents/Simulations/mazia_distr/sim_GEOtop_2011/Q-1.txt');
    
    %  version GEOtop1_1
    % f_tdays	=	1	%	day
    % f_JDfrom0       =		2	%	day
    % f_JD            =		3	%	day
    % f_Q_tot	=	4	%	[m3/s]
    % f_Vsup	=	5	%[m3]
    % f_Vsub	=   6	%[m3]
    % f_Qoutland	=   7   %[m3/s]
    
    %  version GEOtop1_2 and 2_0
    
    % 1)t[days],2)JDfrom0,3)JD,4)Q_tot[m3/s],5)Vsup/Dt[m3/s],6)Vsub/Dt[m3/s],
    % 7)Vchannel[m3], 8)Qoutlandsup[m3/s], 9)Qoutlandsub[m3/s], 10)Qoutbottom[m3/s]
    
    
    f_tdays	=	1;	%	day
    f_JDfrom0       =		2;	%	day
    f_JD            =		3;	%	day
    f_Q_tot	=	4;	%	[m3/s]
    f_Vsup	=	5;	%[m3/s]
    f_Vsub	=   6;	%[m3/s]
    f_channel	=   7;   %[m3]
    f_Qoutland_sup	=   8;  %[m3/s]
    f_Qoutland_sub	=	9;	%	m3/s
    f_Qoutland_bottom	=	10;	%	m3/s
    
end

%% open basin file in basin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(ifbasin)
    disp 'load basin output in basin'
    
    [basin, basin_t, basin_date, basin_labels]=import_GEOtop_2_0_data(basin_name);
    
    if(version~=1.225)
        basin_t=basin(:,2);
    end

    if(version==2.0)
        basin_t=basin(:,1);
    end   
    
    % %  version GEOtop1_32
    % b_JDfrom0       =	1	%
    % b_JD            =	2	%
    % b_t_i           =	3	%	s
    % b_t_f           =	4	%	s
    % b_Prain_bel_can	=	5	%	mm
    % b_Psnow_bel_can	=	6	%	mm
    % b_Prain_abo_can	=	7	%	mm
    % b_Psnow_abo_can	=	8	%	mm
    % b_Tair          =	9	%	C
    % b_Tsurf         =	10	%	C
    % b_Tveg          =	11	%	C
    % b_Evap_surf     =	12	%	mm
    % b_Transp_can	=	13	%	mm
    % b_LE            =	14	%	W/m2
    % b_H             =	15	%	W/m2
    % b_SW            =	16	%	W/m2
    % b_LW            =	17	%	W/m2
    % b_LEv           =	18	%	W/m2
    % b_Hv            =	19	%	W/m2
    % b_SWv           =	20	%	W/m2
    % b_LWv           =	21	%	W/m2
    % b_Swin          =	22	%	W/m2
    % b_Lwin          =	23	%	W/m2
    % b_Mass_bal_err	=	24	%	mm
    
    %  version GEOtop1_1
    % b_tdays	=	1	%	day
    % b_JDfrom0       =		2	%	day
    % b_JD            =		3	%	day
    % b_Prain_bel_can	=	4	%	mm
    % b_Psnow_bel_can	=	5	%	mm
    % b_Prain_abo_can	=	6	%	mm
    % b_Psnow_abo_can	=	7	%	mm
    % b_Tair          =		8	%	C
    % b_Tsurf         =		9	%	C
    % b_Tveg          =		10	%	C
    % b_Evap_surf     =		11	%	mm
    % b_Transp_can	=	12	%	mm
    % b_LE            =		13	%	W/m2
    % b_H             =		14	%	W/m2
    % b_SW            =		15	%	W/m2
    % b_LW            =		16	%	W/m2
    % b_LEv           =		17	%	W/m2
    % b_Hv            =		18	%	W/m2
    % b_SWv           =		19	%	W/m2
    % b_LWv           =		20	%	W/m2
    % b_Swin          =		21	%	W/m2
    % b_Lwin          =		22	%	W/m2
    % b_Mass_bal_err	=	23	%	mm
    
    % Date12[DDMMYYYYhhmm],
    % JulianDayFromYear0[days],TimeFromStart[days],Simulation_Period,Run,
    % 5) Prain_below_canopy[mm],Psnow_below_canopy[mm],Prain_above_canopy[mm],Prain_above_canopy[mm],
    % 9) Tair[C],Tsurface[C],Tvegetation[C],
    % 12) Evap_surface[mm],Transpiration_canopy[mm],LE[W/m2],H[W/m2],SW[W/m2],LW[W/m2],LEv[W/m2],Hv[W/m2],SWv[W/m2],LWv[W/m2],SWin[W/m2],LWin[W/m2],Mass_balance_error[mm]
    
    
    %  version GEOtop1_2		0	%		Date12[DDMMYYYYhhmm]
    b_JDfrom0	=	1	%	day	JulianDayFromYear0[days]
    b_tdays	=	2	%	day	TimeFromStart[days]
    p1_simP	=	3	%	day	Simulation_Period
    p1_run	=	4	%	-	Run
    unknown	=	5	%	-	??
    b_Prain_bel_can	=	6	%	mm	Prain_below_canopy[mm]
    b_Psnow_bel_can	=	7	%	mm	Psnow_below_canopy[mm]
    b_Prain_abo_can	=	8	%	mm	Prain_above_canopy[mm]
    b_Psnow_abo_can	=	9	%	mm	Prain_above_canopy[mm]
    b_Tair          =		10	%	C	Tair[C]
    b_Tsurf         =		11	%	C	Tsurface[C]
    b_Tveg          =		12	%	C	Tvegetation[C]
    b_Evap_surf     =		13	%	mm	Evap_surface[mm]
    b_Transp_can	=	14	%	mm	Transpiration_canopy[mm]
    b_LE            =		15	%	W/m2	LE[W/m2]
    b_H             =		16	%	W/m2	H[W/m2]
    b_SW            =		17	%	W/m2	SW[W/m2]
    b_LW            =		18	%	W/m2	LW[W/m2]
    b_LEv           =		19	%	W/m2	LEv[W/m2]
    b_Hv            =		20	%	W/m2	Hv[W/m2]
    b_SWv           =		21	%	W/m2	SWv[W/m2]
    b_LWv           =		22	%	W/m2	LWv[W/m2]
    b_Swin          =		23	%	W/m2	SWin[W/m2]
    b_Lwin          =		24	%	W/m2	LWin[W/m2]
    b_Mass_bal_err	=	25	%	mm	Mass_balance_error[mm]
    
    if(version==1.225)
        %  version GEOtop1_225		0	%		Date12[DDMMYYYYhhmm]
        b_JDfrom0	=	1	%	day	JulianDayFromYear0[days]
        b_tdays	=	2	%	day	TimeFromStart[days]
        b1_simP	=	3	%	day	Simulation_Period
        b1_run	=	4	%	-	Run
        b_Prain_bel_can	=	5	%	-	Prain_below_canopy[mm]
        b_Psnow_bel_can	=	6	%	mm	Psnow_below_canopy[mm]
        b_Prain_abo_can	=	7	%	mm	Prain_above_canopy[mm]
        b_Psnow_abo_can	=	8	%	mm	Prain_above_canopy[mm]
        b_Pnet	=	9	%	mm	Pnet[mm]
        b_Tair =          		10	%	C	Tair[C]
        b_Tsurf         =		11	%	C	Tsurface[C]
        b_Tveg          =	12	%	C	Tvegetation[C]
        b_Evap_surf     =		13	%	mm	Evap_surface[mm]
        b_Transp_can	=	14	%	mm	Transpiration_canopy[mm]
        b_LE            =		15	%	W/m2	LE[W/m2]
        b_H             =		16	%	W/m2	H[W/m2]
        b_SW            =		17	%	W/m2	SW[W/m2]
        b_LW            =		18	%	W/m2	LW[W/m2]
        b_LEv           =		19	%	W/m2	LEv[W/m2]
        b_Hv            =		20	%	W/m2	Hv[W/m2]
        b_SWv           =		21	%	W/m2	SWv[W/m2]
        b_LWv           =		22	%	W/m2	LWv[W/m2]
        b_Swin          =		23	%	W/m2	SWin[W/m2]
        b_Lwin          =		24	%	W/m2	LWin[W/m2]
        b_Mass_bal_err	=	25	%	mm	Mass_balance_error[mm]
        b_M_Time	=	26	%	s	Mean_Time_Step[s]
    end
    
        
    if(version==2.0)
         %  version GEOtop2_0veg		0	%		Date12[DDMMYYYYhhmm]
        b_JDfrom0	=	1	%	day	JulianDayFromYear0[days]
        b_tdays	=	2	%	day	TimeFromStart[days]
        b1_simP	=	3	%	day	Simulation_Period
        b1_run	=	4	%	-	Run
        b_Prain_bel_can	=	5	%	-	Prain_below_canopy[mm]
        b_Psnow_bel_can	=	6	%	mm	Psnow_below_canopy[mm]
        b_Prain_abo_can	=	7	%	mm	Prain_above_canopy[mm]
        b_Psnow_abo_can	=	8	%	mm	Prain_above_canopy[mm]
        b_Pnet	=	9	%	mm	Pnet[mm]
        b_Tair =          		10	%	C	Tair[C]
        b_Tsurf         =		11	%	C	Tsurface[C]
        b_Tveg          =	12	%	C	Tvegetation[C]
        b_Evap_surf     =		13	%	mm	Evap_surface[mm]
        b_Transp_can	=	14	%	mm	Transpiration_canopy[mm]
        b_LE            =		15	%	W/m2	LE[W/m2]
        b_H             =		16	%	W/m2	H[W/m2]
        b_SW            =		17	%	W/m2	SW[W/m2]
        b_LW            =		18	%	W/m2	LW[W/m2]
        b_LEv           =		19	%	W/m2	LEv[W/m2]
        b_Hv            =		20	%	W/m2	Hv[W/m2]
        b_SWv           =		21	%	W/m2	SWv[W/m2]
        b_LWv           =		22	%	W/m2	LWv[W/m2]
        b_Swin          =		23	%	W/m2	SWin[W/m2]
        b_Lwin          =		24	%	W/m2	LWin[W/m2]
        b_Mass_bal_err	=	25	%	mm	Mass_balance_error[mm]
        b_M_Time	=	26	%	s	Mean_Time_Step[s]
        b_SWE       =   27   %  snow_water_equivalent[mm]
        
    end
    
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load point data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'load multiple point data'


for k=1:npoints
    
    if k<10
        NN=['0',num2str(k)];
    else
        NN=num2str(k);
    end
    %     if (ifpoint)
    %         point_name=[point_root,'.txt'];
    %         Tz_name=[Tz_root,'.txt'];
    %         psiz_name=[psiz_root,'.txt'];
    %         thetaz_name=[thetaz_root,'.txt'];
    %         thetaicez_name=[ thetaicez_root,'.txt'];
    %         snowz_name=[snowz_root,'.txt'];
    %     else
    point_name=[point_root,NN,'.txt'];
    if(ifTz) Tz_name=[Tz_root,NN,'.txt']; end
    if(ifpsiz) psiz_name=[psiz_root,NN,'.txt']; end
    if(ifthetaz) thetaz_name=[thetaz_root,NN,'.txt']; end
    if(ifthetaicez) thetaicez_name=[ thetaicez_root,NN,'.txt']; end
    if(ifsnow) snowz_name=[snowz_root,NN,'.txt']; end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(sprintf('%s in point',point_name))
    
    % GEOtop 1_1 version
    % p1_tdays	=	1	%	day
    % p1_JDfrom0      =		2	%	day
    % p1_JD           =		3	%	day
    % p1_Psnow_bel_can=		4	%	mm
    % p1_Prain_bel_can=		5	%	mm
    %
    % p1_Psnow_abo_can=		6	%	mm
    % p1_Prain_abo_can=		7	%	mm
    % p1_Prain_on_snow=		8	%	mm
    % p1_Wind_s       =		9	%	m/s
    % p1_Wind_d       =		10	%	deg
    % p1_Rh           =		11	%	-
    % p1_Pressure     =		12	%	kPa
    % p1_Tair         =		13	%	C
    % p1_Tdew         =		14	%	C
    % p1_Tsurf        =		15	%	C
    % p1_Tveg         =		16	%	C
    % p1_Tcanopyair =		17	%	C
    % p1_EB_surf      =		18	%	[W/m2]
    % p1_Soil_heat_fl =		19	%	[W/m2]
    % p1_Swin         =		20	%	[W/m2]
    % p1_Swbeam       =		21	%	[W/m2]
    % p1_Swdiff       =		22	%	[W/m2]
    % p1_Lwin         =		23	%	[W/m2]
    % p1_LWin_min     =		24	%	[W/m2]
    % p1_LWin_max     =		25	%	[W/m2]
    % p1_Swnet        =		26	%	[W/m2]
    % p1_Lwnet        =		27	%	[W/m2]
    % p1_H            =		28	%	[W/m2]
    % p1_LE           =		29	%	[W/m2]
    % p1_fc           =		30	%	-
    % p1_LAI          =		31	%	[W/m2]
    % p1_z0veg        =		32	%	m
    % p1_d0veg        =		33	%	m
    % p1_Estored_can =		34	%	[W/m2]
    % p1_SWv          =		35	%	[W/m2]
    % p1_LWv          =		36	%	[W/m2]
    % p1_Hv           =		37	%	[W/m2]
    % p1_LEv          =		38	%	[W/m2]
    % p1_Hg_unveg     =		39	%	[W/m2]
    % p1_LEg_unveg =		40	%	[W/m2]
    % p1_Hg_veg       =		41	%	[W/m2]
    % p1_LEg_veg      =		42	%	[W/m2]
    % p1_Evap_surf	=	43	%	mm
    % p1_Trasp_can	=	44	%	mm
    % p1_Water_on_can	=	45	%	mm
    % p1_Snow_on_can	=	46	%	mm
    % p1_Qveg         =		47	%	-
    % p1_Qsurf        =		48	%	-
    % p1_Qair         =		49	%	-
    % p1_Qcanopyair	=	50	%	-
    % p1_Lobukhov     =		51	%	-
    % p1_Lobukhovcan	=	52	%	m
    % p1_Wind_s_cantop=		53	%	m/s
    % p1_K_deacy_can	=	54	%	-
    % p1_Swup         =		55	%	[W/m2]
    % p1_Lwup         =		56	%	[W/m2]
    % p1_Hup          =		57	%	[W/m2]
    % p1_LEup         =		58	%	[W/m2]
    % p1_Snow_melt	=	59	%	mm
    % p1_Snow_evap	=	60	%	mm
    % p1_Snow_subl	=	61	%	mm
    % p1_Glac_melt	=	62	%	mm
    % p1_Glac_evap	=	63	%	mm
    % p1_Glac_subl	=	64	%	mm
    % p1_depth_thaw	=	65	%	mm
    % p1_depth_wat	=	66	%	mm
    
    % GEOtop 1_2 version and 2_0
    
    p1_JDfrom0      =		1;	%	day
    p1_tdays	=	2;	%	day
    p1_simP = 3;
    p1_run = 4;
    p1_IDpoint = 5;
    p1_Psnow_abo_can=		6;	%	mm
    p1_Prain_abo_can=		7;	%	mm
    p1_Psnow_bel_can=		8;	%	mm
    p1_Prain_bel_can=		9;	%	mm
    p1_Prain_on_snow=		10;	%	mm
    p1_Wind_s       =		11;	%	m/s
    p1_Wind_d       =		12;	%	deg
    p1_Rh           =		13;	%	-
    p1_Pressure     =		14;	%	kPa
    p1_Tair         =		15;	%	C
    p1_Tdew         =		16;	%	C
    p1_Tsurf        =		17;	%	C
    p1_Tveg         =		18;	%	C
    p1_Tcanopyair =		19;	%	C
    p1_EB_surf      =		20;	%	[W/m2]
    p1_Soil_heat_fl =		21;	%	[W/m2]
    p1_Swin         =		22;	%	[W/m2]
    p1_Swbeam       =		23;	%	[W/m2]
    p1_Swdiff       =		24;	%	[W/m2]
    p1_Lwin         =		25;	%	[W/m2]
    p1_LWin_min     =		26;	%	[W/m2]
    p1_LWin_max     =		27;	%	[W/m2]
    p1_Swnet        =		28;	%	[W/m2]
    p1_Lwnet        =		29;	%	[W/m2]
    p1_H            =		30;	%	[W/m2]
    p1_LE           =		31;	%	[W/m2]
    p1_fc           =		32;	%	-
    p1_LAI          =		33;	%	[W/m2]
    p1_z0veg        =		34;	%	m
    p1_d0veg        =		35;	%	m
    p1_Estored_can =		36;	%	[W/m2]
    p1_SWv          =		37;	%	[W/m2]
    p1_LWv          =		38;	%	[W/m2]
    p1_Hv           =		39;	%	[W/m2]
    p1_LEv          =		40;	%	[W/m2]
    p1_Hg_unveg     =		41;	%	[W/m2]
    p1_LEg_unveg =		42;	%	[W/m2]
    p1_Hg_veg       =		43;	%	[W/m2]
    p1_LEg_veg      =		44;	%	[W/m2]
    p1_Evap_surf	=	45;	%	mm
    p1_Trasp_can	=	46;	%	mm
    p1_Water_on_can	=	47;	%	mm
    p1_Snow_on_can	=	48;	%	mm
    p1_Qveg         =		49;	%	-
    p1_Qsurf        =		50;	%	-
    p1_Qair         =		51;	%	-
    p1_Qcanopyair	=	52;	%	-
    p1_Lobukhov     =		53;	%	-
    p1_Lobukhovcan	=	54;	%	m
    p1_Wind_s_cantop=		55;	%	m/s
    p1_K_decay_can	=	56;	%	-
    p1_Swup         =		57;	%	[W/m2]
    p1_Lwup         =		58;	%	[W/m2]
    p1_Hup          =		59;	%	[W/m2]
    p1_LEup         =		60;	%	[W/m2]
    p1_Snow_depth = 61;  %   mm
    p1_Snow_swe = 62;  %   mm
    p1_Snow_density = 63;  %   mm
    p1_Snow_temp = 64; %   mm
    p1_Snow_melt	=	65;	%	mm
    p1_Snow_subl	=	66;	%	mm
    p1_Snow_blown	=	67;	%	mm
    p1_Snow_sub_blown	=	68;	%	mm
    p1_Glac_depth	=	69;	%	mm
    p1_Glac_swe	=	70;	%	mm
    p1_Glac_density	=	71;	%	mm
    p1_Glac_temperature	=	72;	%	mm
    p1_Glac_melt	=	73;	%	mm
    p1_Glac_subl	=	74;	%	mm
    p1_depth_thaw	=	75;	%	mm
    p1_depth_wat	=	76;	%	mm
    
    if(or((version==1.225),(version==2.0)))
        % GEOtop 1_225 and 2_0 version			%		0	Date12[DDMMYYYYhhmm]
        p1_low_thaw	=	75	%	mm	75	lowest_thawed_soil_depth[mm]
        p1_high_thaw	=	76	%	mm	76	highest_thawed_soil_depth[mm]
        p1_low_wat	=	77	%	mm	77	lowest_water_table_depth[mm]
        p1_high_wat	=	78	%	mm	78	highest_water_table_depth[mm]
        
    end
    
    
    [tmp, point_t, point_date, point_labels]=import_GEOtop_2_0_data(point_name);
    if(k==1), point=zeros(size(tmp,1),size(tmp,2),npoints); end;
    point(:,:,k)=tmp;
    clear tmp
    
    point( point==-9999 )=NaN;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Tz file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(ifTz)
        disp(sprintf('%s in Tz',Tz_name))
        % Date12[DDMMY6YYhhmm], 1) JulianDayFromYear0[days], 2) TimeFromStart[days], 3) Simulation_Period,Run, 4) IDpoint,
        % 5)'Layer1 ' 4)'Layer2 ' 7)'... 'Layer Nz'
        [tmp, Tz_t, Tz_date, Tz_labels]=import_GEOtop_2_0_data(Tz_name);
        if(k==1), Tz=zeros(size(tmp,1),size(tmp,2),npoints); end;
        Tz(:,:,k)=tmp;
        clear tmp
        Tz( Tz==-9999 )=NaN;
        
        soiltem_file= importdata(Tz_name, ',',1);
        header=soiltem_file.textdata(1,7:end);
        header_num=str2num(cell2mat(header))/1000; %in m
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % psiz file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(ifpsiz)
        disp(sprintf('%s in psiz',psiz_name))
        % Date12[DDMMY6YYhhmm], 1) JulianDayFromYear0[days], 2) TimeFromStart[days], 3) Simulation_Period,Run, 4) IDpoint,
        % 5)'Layer1 ' 4)'Layer2 ' 7)'... 'Layer Nz'
        [tmp, psiz_t, psiz_date, psiz_labels]=import_GEOtop_2_0_data(psiz_name);
        if(k==1), psiz=zeros(size(tmp,1),size(tmp,2),npoints); end;
        psiz(:,:,k)=tmp; clear tmp
        psiz( psiz==-9999 )=NaN;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % thetaz file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (ifthetaz)
        disp(sprintf('%s in thetaz',thetaz_name))
        % Date12[DDMMY6YYhhmm], 1) JulianDayFromYear0[days], 2) TimeFromStart[days], 3) Simulation_Period,Run, 4) IDpoint,
        % 5)'Layer1 ' 4)'Layer2 ' 7)'... 'Layer Nz'
        [tmp, thetaz_t, thetaz_date, thetaz_labels]=import_GEOtop_2_0_data(thetaz_name);
        if(k==1), thetaz=zeros(size(tmp,1),size(tmp,2),npoints); end;
        thetaz(:,:,k)=tmp; clear tmp
        thetaz( thetaz==-9999 )=NaN;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % thetaicez file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (ifthetaicez)
        disp(sprintf('%s in thetaicez',thetaicez_name))
        % Date12[DDMMY6YYhhmm], 1) JulianDayFromYear0[days], 2) TimeFromStart[days], 3) Simulation_Period,Run, 4) IDpoint,
        % 5)'Layer1 ' 4)'Layer2 ' 7)'... 'Layer Nz'
        [tmp, thetaicez_t, thetaicez_date, thetaicez_labels]=import_GEOtop_2_0_data(thetaicez_name);
        if(k==1), thetaicez=zeros(size(tmp,1),size(tmp,2),npoints); end;
        thetaicez(:,:,k)=tmp; clear tmp
        thetaicez( thetaicez==-9999 )=NaN;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % snowz file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BUG : still to update to 2.0 version
    if (ifsnow)
        disp(sprintf('%s in snowz',snowz_name))
        
        % % GEOtop 1_1 version
        % s1_tdays	=	1	%	day
        % s1_JDfrom0	=	2	%	day
        % s1_JD	=	3	%	day
        % s1_Snowtype	=	4	%	-
        % s1_Snowdepth	=	5	%	mm
        % s1_SWE	=	6	%	mm
        % s1_Temp_aver	=	7	%	C
        % s1_Density_aver	=	8	%	kg/m3
        % s1_melting	=	9	%	mm
        % s1_sublimation	=	10	%	mm
        % s1_evaporation	=	11	%	mm
        % s1_nlayer	=	12	%	-
        % s1_Bstrans	=	13	%	mm
        % s1_Bssubl	=	14	%	mm
        % s1_Bstot	=	15	%	mm
        % s1_Dlayer	=	[16:15+Ns]	%	mm
        % s1_Zlayer	=	[15+Ns+1:15+Ns*2]	%	mm
        % s1_Temp	=	[15+Ns*2+1:15+Ns*3]	%	C
        % s1_Density	=	[15+Ns*3+1:15+Ns*4]	%	[kg/m3]
        % s1_SWE_l	=	[15+Ns*4+1:15+Ns*5]	%	mm
        % s1_ice	=	[15+Ns*5+1:15+Ns*6]	%	[kg/m2]
        % s1_water	=	[15+Ns*6+1:15+Ns*7]	%	[kg/m2
        
        % GEOtop 1_2 version
        s1_JDfrom0	=	1	%	day
        s1_tdays	=	2	%	day
        s1_simP	=	3	%	-
        s1_run	=	4	%	-
        s1_IDpoint	=	5	%	-
        s1_SWE_l	=	[6:5+Ns]	%	mm
        s1_Snowdepth	=	[5+Ns+1:5+Ns*2]	%	mm
        s1_Density	=	[5+Ns*2+1:5+Ns*3]	%	kg/m3
        s1_Temp	=	[5+Ns*3+1:5+Ns*4]	%	C
        s1_ice	=	[5+Ns*4+1:5+Ns*5]	%	-
        s1_water	=	[5+Ns*5+1:5+Ns*6]	%	-
        
        [tmp, snowz_t, snowz_date, snowz_labels]=import_GEOtop_1_2_data(snowz_name);
        if(k==1),  snowz=zeros(size(tmp,1),size(tmp,2),npoints); end;
        snowz(:,:,k)=tmp; clear tmp
        snowz( snowz==-9999 )=NaN;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calcolates useful quantities
dts=Dt_out_pix*3600; % dt ouptut in s
dth=Dt_out_pix; % dt ouptut in h
dtd=Dt_out_pix/24; % dt ouptut in d

% calculate LEup ETup for version1.225-9
% Hup=(H+fc*Hv)
% H= fc * Hgveg + (1-fc)* Hgunveg 
% Hup = fc * (Hgveg + Hv ) + (1-fc) * Hgunveg 
% p1_H p1_fc p1_Hv p1_Hg_unveg p1_Hg_veg
%
% LEup= (LE+fc*LEv)
% LE =fc * LEgveg + (1-fc)* LEgunveg 
% LEup = fc * (LEgveg + LEv ) + (1-fc) * LEgunveg 
% p1_LE p1_fc p1_LEv p1_LEg_unveg p1_LEg_veg

% For those versions total LE and H above the canopy are calculated in postprocess
if(or((version==1.225),(version==2.0)))
    
np=size(point,1);
rqv=zeros(np,npoints);
rqg=zeros(np,npoints);
rqa=zeros(np,npoints);
LEup_veg=zeros(np,npoints);
LEup=zeros(np,npoints);
    for k=1:npoints
%         point(:,p1_Hup,k)= point(:,p1_H,k)+point(:,p1_fc,k).*point(:,p1_Hv,k);
%         point(:,p1_LEup,k)= point(:,p1_LE,k)+point(:,p1_fc,k).*point(:,p1_LEv,k);
        point(:,p1_Hup,k)= point(:,p1_fc,k).*(point(:,p1_Hg_veg,k)+point(:,p1_Hv,k))+...
            (1-point(:,p1_fc,k)).*point(:,p1_Hg_unveg,k);
        point(:,p1_LEup,k)= point(:,p1_fc,k).*(point(:,p1_LEg_veg,k)+point(:,p1_LEv,k))+...
            (1-point(:,p1_fc,k)).*point(:,p1_LEg_unveg,k);
   
        
% 2015.06.04: updated way to calculate LEup and Hup
%         Recalculate LE_up using the parallel resistance approach
%         LEv=(qv-qca)/rv   ->rv
%         LEgv=(qs-qca)/rg   -> rg
%         ra=(rv*rg)/(rv+rg) -> ra
%         LEup_veg=(qca-qa)/ra -> LEup_veg
%         LEup=fc*LEup_veg+(1-fc)*LEg

%   Names for LE
%   LEv p1_LEv
%   qv p1_Qveg
%   qca p1_Qcanopyair
%   qs p1_Qsurf
%   qa p1_Qair
% p1_LEg_veg
% p1_LEg_unveg
% p1_fc
% p1_LEup



rqv(:,k)=(point(:,p1_Qveg,k)-point(:,p1_Qcanopyair,k))./point(:,p1_LEv,k);
rqg(:,k)=(point(:,p1_Qsurf,k)-point(:,p1_Qcanopyair,k))./point(:,p1_LEg_veg,k);
rqa(:,k)=rqv(:,k).*rqg(:,k)./(rqv(:,k)+rqg(:,k));
LEup_veg(:,k)=(point(:,p1_Qcanopyair,k)-point(:,p1_Qair,k))./rqa(:,k);
LEup(:,k)=point(:,p1_fc,k).*LEup_veg(:,k)+(1-point(:,p1_fc,k)).*point(:,p1_LEg_unveg,k);

%   Names for H
%   Hv p1_Hv
%   Ta p1_Tair     
%   Ts p1_Tsurf        
%   Tv p1_Tveg         
%   Tca p1_Tcanopyair 
% p1_Hg_veg
% p1_Hg_unveg
% p1_fc
% p1_Hup      
        
  

   end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load soil data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load _soil.txt file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Each block corresponds to each soil type given in the correspondent map.
% If no soil type file is provided, the first block parameters and initial conditions will be used
% The following data give the values of the main soil properties (column) to each soil layer (row).

% 1)Dz,2)Kh,3)Kv,4)vwc_r,5)vwc_w,6)vwc_fc,7)vwc_s,8)alpha,9)n,10)stor
%
% layer 1
% layer 2
% layer 3
% .....
% #1	Thickness(mm)
% #2	Kh(mm/s)
% #3	Kv(mm/s)
% #4	residual water content
% #5  water content at the wilting point
% #6  field capacity water content
% #7  saturated water content
% #8	alpha(mm^-1)
% #9	n
% #10	Specific storage (mm^-1)

if(ifsoil)
    disp 'load _soil.txt in soilpar'
    
    for k=1:n_soils
        if k<10
            NN=['0',num2str(k)];
        else
            NN=num2str(k);
        end
        soil_name=[soil_root,NN,'.txt'];
        
        disp(sprintf('%s in soil',soil_name))
        soilpar= importdata(soil_name, ',',1);
        
        % 1)Dz, 2)Kh, 3)Kv, 4)vwc_r, 5)vwc_w, 6)vwc_fc, 7)vwc_s, 8)alpha, 9)n, 10)stor
        theta_r(:,k)=soilpar.data(:,4);
        theta_s(:,k)=soilpar.data(:,7);
    end
    Nz=size(soilpar.data,1); % number of vertical layers
    Dz=soilpar.data(:,1);
    
else
    % disp 'load _soil.txt in soilpar'
    disp 'Soil parameters from geotop.inpts'
    
    searchterm='SoilLayerThicknesses' ;
    for i=1:n
        if strcmp(wordc{i},searchterm)==1;
            tmp1= wordc{i+2};
            tmp2=strrep(tmp1, ',', ';');
            Dz=str2num(tmp2);
            clear 'tmp1' 'tmp2'
        end
    end
    Nz=size(Dz,1);
    
    searchterm='ThetaSat' ;
    for i=1:n
        if strcmp(wordc{i},searchterm)==1;
            tmp1= wordc{i+2};
            tmp2=strrep(tmp1, ',', ';');
            theta_s0=str2num(tmp2);
            clear 'tmp1' 'tmp2'
        end
    end
    
    if(size(theta_s0,1)<Nz);
        theta_s=ones(Nz,1)*theta_s0;
    else
        theta_s=theta_s0;
    end
        
    searchterm='ThetaRes' ;
    for i=1:n
        if strcmp(wordc{i},searchterm)==1;
            tmp1= wordc{i+2};
            tmp2=strrep(tmp1, ',', ';');
            theta_r0=str2num(tmp2);
            clear 'tmp1' 'tmp2'
        end
    end
    
    if(size(theta_r0,1)<Nz);
        theta_r=ones(Nz,1)*theta_r0;
    else
        theta_r=theta_r0;
    end
    
end

% missed part with soil parameters in geotop.inpts
%     if (ifpoint)
%         soil_name=inpts_name;
%         % '1)SoilLayerThicknesses'
%         searchterm='SoilLayerThicknesses' ;
%         for i=1:n
%             if strcmp(wordc{i},searchterm)==1;
%                 index=i;
% %                 SoilLayerThicknesses= str2double(wordc{i+2})
%             end
%         end
%
%         j=index
%         j_1=1;
%         flag=0;
%         while (flag==0)
%             if (~isnan(str2double(wordc{j+2})))
%                 Dz(j_1)= str2double(wordc{j+2});
%                 j=j+1;
%                 j_1=j_1+1;
%             else
%                 flag=1
%             end
%         end
%
%         Nz=size(Dz,2); % number of vertical layers
%         % '1)ThetaRes'
%         searchterm='ThetaRes' ;
%         for i=1:n
%             if strcmp(wordc{i},searchterm)==1;
%                 index=i;
%             end
%         end
%
%         j=index
%         j_1=1;
%         flag=0;
%         while (flag==0)
%             if (~isnan(str2double(wordc{j+2})))
%                 theta_r(j_1)= str2double(wordc{j+2});
%                 j=j+1;
%                 j_1=j_1+1;
%             else
%                 flag=1
%             end
%         end
%         % '1)ThetaRes'
%         searchterm='ThetaSat' ;
%         for i=1:n
%             if strcmp(wordc{i},searchterm)==1;
%                 index=i;
%             end
%         end
%
%         j=index
%         j_1=1;
%         flag=0;
%         while (flag==0)
%             if (~isnan(str2double(wordc{j+2})))
%                 theta_s(j_1)= str2double(wordc{j+2});
%                 j=j+1;
%                 j_1=j_1+1;
%             else
%                 flag=1
%             end
%         end
%     else
%

%% calculates saturation
% bug: theta_r(k,h) the index h shold be the soil type of the outptut point
% h
if (ifthetaz)
    if(ifsatz)
    satz=thetaz;
    for h=1:npoints
    for k=1:Nz
        satz(:,k+noffset,h)=(thetaz(:,k+noffset,h)-theta_r(k,h))./(theta_s(k,h)-theta_r(k,h));
    end
    end
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unload simulation
unload_GEOtop_2_0_sim_path(mysim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculates point  energy budget components
%Sw_net_atm=p1_Swin - p1_Swup
Sw_net_atm=point(:,p1_Swin) - point(:,p1_Swup);
%Lw_net_atm=p1_Lwin - p1_Lwup
Lw_net_atm=point(:,p1_Lwin) - point(:,p1_Lwup);
Rn_net_atm=Sw_net_atm+Lw_net_atm;
%p1_Hup
%p1_LEup
%EB_ab_veg=Rn_net_atm-p1_Hup-p1_LEup
EB_ab_veg=Rn_net_atm-point(:,p1_Hup)-point(:,p1_LEup);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% averages results over a daily basis: flow basin point Tz thetaz psiz
if (daily)
    disp 'averages results over a daily basis'
    
    JDfrom0_p1=point(:,p1_JDfrom0,1);
    
    if (~ifpoint)
        JDfrom0_f = flows(:,f_JDfrom0);
    end
    
    JDfrom0_s1 = snowz(:,s1_JDfrom0);
    % JDfrom0_t1=point(:,t1_JDfrom0,1);
    
    MNTH_LNGTH_LEAP=[31 29 31 30 31 30 31 31 30 31 30 31];
    MNTH_LNGTH=[31 29 31 30 31 30 31 31 30 31 30 31];
    date1=datevec(JDfrom0_p1(1));
    Y1=date1(1);
    date2=datevec(JDfrom0_p1(end));
    Y2=date2(1);
    
    if (ifpoint == 0)
        flow_ts = timeseries(flows,JDfrom0_f);
    end
    
    point_ts = timeseries(point,JDfrom0_p1);
    
    if (ifTz)
        Tz_ts = timeseries(Tz(:,1+noffset:end),Tz(:,1));
    end
    
    if(ifpsiz)
        psiz_ts = timeseries(psiz(:,1+noffset:end),psiz(:,1));
    end
    
    if(ifpsiztot)
        psiztot_ts = timeseries(psiztot(:,1+noffset:end),psiztot(:,1));
    end
    
    if (ifthetaz)
        thetaz_ts = timeseries(thetaz(:,1+noffset:end),thetaz(:,1));
    end
    
    if (ifthetaicez)
        thetaicez_ts = timeseries(thetaicez(:,1+noffset:end),thetaicez(:,1));
    end
    
    snowz_ts = timeseries(snowz,JDfrom0_s1);
    
    Sw_net_atm_ts = timeseries(Sw_net_atm,JDfrom0_p1);
    
    Lw_net_atm_ts = timeseries(Lw_net_atm,JDfrom0_p1);
    
    Rn_net_atm_ts = timeseries(Rn_net_atm,JDfrom0_p1);
    
    EB_ab_veg_ts = timeseries(EB_ab_veg,JDfrom0_p1);
    
    
    % Calculate daily means
    i=1;
    for YEAR=Y1:Y2
        for MONTH=1:12
            % if leap year
            if ((YEAR==1980) || (YEAR==1984) || (YEAR==1988) || (YEAR==1992) || (YEAR==1996) || (YEAR==2000) || (YEAR==2004) || (YEAR==2008))
                for DAY = 1:MNTH_LNGTH_LEAP(MONTH)
                    if ((datenum(YEAR,MONTH,DAY)>=datenum(date1)) && (datenum(YEAR,MONTH,DAY)<=datenum(date2)))
                        starttime=datenum(YEAR,MONTH,DAY,0,0,0);
                        endtime=datenum(YEAR,MONTH,DAY,23,59,59);
                        if (ifpoint == 0)
                            flow_daily_ts = getsampleusingtime(flow_ts, starttime, endtime);
                            flow_d(i,:)=mean(flow_daily_ts.data,1);
                        end
                        
                        point_daily_ts=getsampleusingtime(point_ts, starttime, endtime);
                        point_d(i,:)=mean(point_daily_ts.data,1);
                        
                        if (ifTz)
                            Tz_daily_ts=getsampleusingtime(Tz_ts, starttime, endtime);
                            Tz_d(i,:)=mean(Tz_daily_ts.data,1);
                        end
                        if(ifpsiz)
                            psiz_daily_ts=getsampleusingtime(psiz_ts, starttime, endtime);
                            psiz_d(i,:)=mean(psiz_daily_ts.data,1);
                        end
                        
                        if(ifpsiztot)
                            psiztot_daily_ts=getsampleusingtime(psiztot_ts, starttime, endtime);
                            psiztot_d(i,:)=mean(psiztot_daily_ts.data,1);
                        end
                        
                        if (ifthetaz)
                            thetaz_daily_ts=getsampleusingtime(thetaz_ts, starttime, endtime);
                            thetaz_d(i,:)=mean(thetaz_daily_ts.data,1);
                        end
                        
                        if (ifthetaicez)
                            thetaicez_daily_ts=getsampleusingtime(thetaicez_ts, starttime, endtime);
                            thetaicez_d(i,:)=mean(thetaicez_daily_ts.data,1);
                        end
                        
                        snowz_daily_ts=getsampleusingtime(snowz_ts, starttime, endtime);
                        snowz_d(i,:)=mean(snowz_daily_ts.data,1);
                        
                        Sw_net_atm_daily_ts=getsampleusingtime(Sw_net_atm_ts, starttime, endtime);
                        Sw_net_atm_d(i,:)=mean(Sw_net_atm_daily_ts.data,1);
                        
                        Lw_net_atm_daily_ts=getsampleusingtime(Lw_net_atm_ts, starttime, endtime);
                        Lw_net_atm_d(i,:)=mean(Lw_net_atm_daily_ts.data,1);
                        
                        Rn_net_atm_daily_ts=getsampleusingtime(Rn_net_atm_ts, starttime, endtime);
                        Rn_net_atm_d(i,:)=mean(Rn_net_atm_daily_ts.data,1);
                        
                        EB_ab_veg_daily_ts=getsampleusingtime(EB_ab_veg_ts, starttime, endtime);
                        EB_ab_veg_d(i,:)=mean(EB_ab_veg_daily_ts.data,1);
                        
                        daily_date(i)=datenum(YEAR,MONTH,DAY,0,0,0);
                        i=i+1;
                    end
                end
            else
                for DAY = 1:MNTH_LNGTH(MONTH)
                    if ((datenum(YEAR,MONTH,DAY)>=datenum(date1)) && (datenum(YEAR,MONTH,DAY)<=datenum(date2)))
                        starttime=datenum(YEAR,MONTH,DAY,0,0,0);
                        endtime=datenum(YEAR,MONTH,DAY,23,59,59);
                        
                        if (ifpoint == 0)
                            flow_daily_ts = getsampleusingtime(flow_ts, starttime, endtime);
                            flow_d(i,:)=mean(flow_daily_ts.data,1);
                        end
                        
                        point_daily_ts=getsampleusingtime(point_ts, starttime, endtime);
                        point_d(i,:)=mean(point_daily_ts.data,1);
                        
                        if (ifTz)
                            Tz_daily_ts=getsampleusingtime(Tz_ts, starttime, endtime);
                            Tz_d(i,:)=mean(Tz_daily_ts.data,1);
                        end
                        if(ifpsiz)
                            psiz_daily_ts=getsampleusingtime(psiz_ts, starttime, endtime);
                            psiz_d(i,:)=mean(psiz_daily_ts.data,1);
                        end
                        
                        if(ifpsiztot)
                            psiztot_daily_ts=getsampleusingtime(psiztot_ts, starttime, endtime);
                            psiztot_d(i,:)=mean(psiztot_daily_ts.data,1);
                        end
                        
                        if (ifthetaz)
                            thetaz_daily_ts=getsampleusingtime(thetaz_ts, starttime, endtime);
                            thetaz_d(i,:)=mean(thetaz_daily_ts.data,1);
                        end
                        
                        if (ifthetaicez)
                            thetaicez_daily_ts=getsampleusingtime(thetaicez_ts, starttime, endtime);
                            thetaicez_d(i,:)=mean(thetaicez_daily_ts.data,1);
                        end
                        
                        snowz_daily_ts=getsampleusingtime(snowz_ts, starttime, endtime);
                        snowz_d(i,:)=mean(snowz_daily_ts.data,1);
                        
                        Sw_net_atm_daily_ts=getsampleusingtime(Sw_net_atm_ts, starttime, endtime);
                        Sw_net_atm_d(i,:)=mean(Sw_net_atm_daily_ts.data,1);
                        
                        Lw_net_atm_daily_ts=getsampleusingtime(Lw_net_atm_ts, starttime, endtime);
                        Lw_net_atm_d(i,:)=mean(Lw_net_atm_daily_ts.data,1);
                        
                        Rn_net_atm_daily_ts=getsampleusingtime(Rn_net_atm_ts, starttime, endtime);
                        Rn_net_atm_d(i,:)=mean(Rn_net_atm_daily_ts.data,1);
                        
                        EB_ab_veg_daily_ts=getsampleusingtime(EB_ab_veg_ts, starttime, endtime);
                        EB_ab_veg_d(i,:)=mean(EB_ab_veg_daily_ts.data,1);
                        
                        daily_date(i)=datenum(YEAR,MONTH,DAY,0,0,0);
                        i=i+1;
                    end
                end
            end
        end
    end
end
% Calculate Annual Mean
%     i=1;
%     for YEAR=Y1:Y2
%         starttime=datenum(YEAR,1,1,0,0,0);
%         endtime=datenum(YEAR,12,31,23,59,59);
%         Tz_annual_ts(i) = getsampleusingtime(Tz_ts, starttime, endtime)
%         Tz_annual_date(i)=datenum(YEAR,6,30,0,0,0);
%         i=i+1;
%     end
%
%     for i=1:length(Tz_annual_ts)
%         Tz_annual_mean(i,:)=mean(Tz_annual_ts(i));
%     end


if(ifbasin)
    %     if(Dt_out_bas<=24)
    %         basin_d=resample(basin,24/Dt_out_bas);
    %     else
    %         clear basin_d
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% calculates basin averaged energy budget components
    b_EB=basin(:,b_SW)+basin(:,b_LW)-basin(:,b_LE)-basin(:,b_H);
    b_E_stor_can=basin(:,b_SWv) + basin(:,b_LWv) - basin(:,b_Hv) - basin(:,b_LEv);
    % SW_b=basin(:,8); LW_b=basin(:,9);
    % Rn_b=SW_b+LW_b;
    % H_b=basin(:,10); LE_b=basin(:,11); G_b=Rn_b-H_b-LE_b;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates total basin storage
% SS=30)SSup[mm]+ 31)SSub[mm]+ 32)wt[mm]+ 33)SWE[mm]
%SS=sum(basin(:,30:33),2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basin cumulative mass budget
disp 'calculates basin cumulative mass budget'
% basin_c=cumsum(basin);
% E_c=sum(basin_c(:,5:7),2); % 5)Esoil[mm] 6)Ecv[mm] 7)Etc[mm]
% P_c=sum(basin_c(:,20:21),2); % 20)Prain[mm], 21)Psnow[mm]
% R_c=basin_c(:,29); % 29)R_tot[mm]
% S_c=P_c-E_c-R_c; % total cumulative difference mass budget

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basin cumulative energy budget
disp 'calculates basin cumulative energy budget'
% Rn_c=cumsum(Rn_b); LE_c=cumsum(LE_b); H_c=cumsum(H_b); G_c=cumsum(G_b);
% Eb_err_c=Rn_c-LE_c-H_c-G_c;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simfile=[mysim,simname,'.mat'];
fprintf('saves output in %s\n',simfile)
save(simfile);







