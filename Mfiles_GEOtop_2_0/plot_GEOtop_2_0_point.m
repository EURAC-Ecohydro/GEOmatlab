%% plot_GEOtop_2_0_point.m
% Copyright, 2010 Giacomo Bertoldi;
%
% plots GEOTOP 2_0 temporal outputs for a chosen control point ;
% load before load_GEOtop_2_0.m
% you need to define the chosen point # in the script
%%
if(ifbasin==1)
	tbeg_bas=min(basin_t); tend_bas=max(basin_t);
end
tbeg_pnt=min(point_t); tend_pnt=max(point_t);

% p_start=datenum('02-Jan-2009 00:00:00');
% p_end=datenum('05-Jan-2009 00:00:00');
% 
% p_start=datenum('02-Jun-2009 00:00:00');
% p_end=datenum('05-Jun-2009 00:00:00');

%% define point to be plotted if you have multiple point ouptuts

% choose outptut point to be plotted
p=1;

if(p>npoints) 
    disp(sprintf('Error: max number of outputpoints possible = %d',npoints))
    return
end

disp(sprintf('Simulation %s',mysim));
disp(sprintf('Ouptut plotted on %s',date));
disp(sprintf('Plots from %s to %s',datestr(tbeg_pnt),datestr(tend_pnt)));
disp(sprintf('Plotted output point number %d',p));

%% 'po03_point0001.txt in point'
disp(sprintf('%s in point',point_name))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% new labels GEOtop 1_2 version
%     p1_JDfrom0      =		1;	%	day
%     p1_tdays	=	2;	%	day
%     p1_simP = 3;
%     p1_run = 4;
%     p1_IDpoint = 5;
%     p1_Psnow_abo_can=		6;	%	mm
%     p1_Prain_abo_can=		7;	%	mm
%     p1_Psnow_bel_can=		8;	%	mm
%     p1_Prain_bel_can=		9;	%	mm
%     p1_Prain_on_snow=		10;	%	mm
%     p1_Wind_s       =		11;	%	m/s
%     p1_Wind_d       =		12;	%	deg
%     p1_Rh           =		13;	%	-
%     p1_Pressure     =		14;	%	kPa
%     p1_Tair         =		15;	%	C
%     p1_Tdew         =		16;	%	C
%     p1_Tsurf        =		17;	%	C
%     p1_Tveg         =		18;	%	C
%     p1_Tcanopyair =		19;	%	C
%     p1_EB_surf      =		20;	%	[W/m2]
%     p1_Soil_heat_fl =		21;	%	[W/m2]
%     p1_Swin         =		22;	%	[W/m2]
%     p1_Swbeam       =		23;	%	[W/m2]
%     p1_Swdiff       =		24;	%	[W/m2]
%     p1_Lwin         =		25;	%	[W/m2]
%     p1_LWin_min     =		26;	%	[W/m2]
%     p1_LWin_max     =		27;	%	[W/m2]
%     p1_Swnet        =		28;	%	[W/m2]
%     p1_Lwnet        =		29;	%	[W/m2]
%     p1_H            =		30;	%	[W/m2]
%     p1_LE           =		31;	%	[W/m2]
%     p1_fc           =		32;	%	-
%     p1_LAI          =		33;	%	[W/m2]
%     p1_z0veg        =		34;	%	m
%     p1_d0veg        =		35;	%	m
%     p1_Estored_can =		36;	%	[W/m2]
%     p1_SWv          =		37;	%	[W/m2]
%     p1_LWv          =		38;	%	[W/m2]
%     p1_Hv           =		39;	%	[W/m2]
%     p1_LEv          =		40;	%	[W/m2]
%     p1_Hg_unveg     =		41;	%	[W/m2]
%     p1_LEg_unveg =		42;	%	[W/m2]
%     p1_Hg_veg       =		43;	%	[W/m2]
%     p1_LEg_veg      =		44;	%	[W/m2]
%     p1_Evap_surf	=	45;	%	mm
%     p1_Trasp_can	=	46;	%	mm
%     p1_Water_on_can	=	47;	%	mm
%     p1_Snow_on_can	=	48;	%	mm
%     p1_Qveg         =		49;	%	-
%     p1_Qsurf        =		50;	%	-
%     p1_Qair         =		51;	%	-
%     p1_Qcanopyair	=	52;	%	-
%     p1_Lobukhov     =		53;	%	-
%     p1_Lobukhovcan	=	54;	%	m
%     p1_Wind_s_cantop=		55;	%	m/s
%     p1_K_deacy_can	=	56;	%	-
%     p1_Swup         =		57;	%	[W/m2]
%     p1_Lwup         =		58;	%	[W/m2]
%     p1_Hup          =		59;	%	[W/m2]
%     p1_LEup         =		60;	%	[W/m2]
%     p1_Snow_depth = 61;  %   mm
%     p1_Snow_swe = 62;  %   mm
%     p1_Snow_density = 63;  %   mm
%     p1_Snow_temp = 64; %   mm
%     p1_Snow_melt	=	65;	%	mm
%     p1_Snow_subl	=	66;	%	mm
%     p1_Snow_blown	=	67;	%	mm
%     p1_Snow_sub_blown	=	68;	%	mm
%     p1_Glac_depth	=	69;	%	mm
%     p1_Glac_swe	=	70;	%	mm
%     p1_Glac_density	=	71;	%	mm
%     p1_Glac_temperature	=	72;	%	mm
%     p1_Glac_melt	=	73;	%	mm
%     p1_Glac_subl	=	74;	%	mm
%     p1_depth_thaw	=	75;	%	mm
%     p1_depth_wat	=	76;	%	mm
 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 11: plot temperatures'

% p1_Tair         =	15	%	C
% p1_Tsurf        =	17	%	C
% p1_Tveg         =	18	%	C
% p1_Tcanopyair	=	19	%	C

% clean point(:,p1_Tcanopyair)
Tcanopyair=point(:,p1_Tcanopyair,p);
tmp=find(Tcanopyair<-50);
Tcanopyair(tmp)=NaN;

figure(11),clf
plot(point_t,point(:,p1_Tair,p)), hold on
plot(point_t,Tcanopyair,'k')
plot(point_t,point(:,p1_Tveg,p),'g')
plot(point_t,point(:,p1_Tsurf,p),'r')
legend(' Tair[C]',' Ts canopy air',' T veg',' T surf ')
title('temperatures [C]','fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks');
dynamicDateTicks

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 111: plot meteo input'

% p1_Tair         =	15	%	C
% p1_Tdew        =	16	%	C
% p1_Wind_s         =	11	%	C
% p1_Rh	=	13	%	C


figure(111),clf
ax(1)=subplot(3,1,1)
plot(point_t,point(:,p1_Tair,p),'r'), hold on
plot(point_t,point(:,p1_Tdew,p),'b')
legend(' Tair[C]',' Tdew[C]')
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks');
dynamicDateTicks
ax(2)=subplot(3,1,2)
plot(point_t,point(:,p1_Rh,p),'b')
legend(' RH[-]')
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks');
dynamicDateTicks
ax(3)=subplot(3,1,3)
plot(point_t,point(:,p1_Wind_s,p),'k')
legend(' Wind[m/s]')
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks');
dynamicDateTicks
linkaxes(ax,'x');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 12: plot SW rad'
%%%
% p1_Swin         =	22	%	[W/m2]
% p1_SWv          =	37	%	[W/m2]
% p1_Swnet        =	28	%	[W/m2]
% p1_Swup         =	58	%	[W/m2]
figure(12),clf
plot(point_t,point(:,p1_Swin,p)), hold on
plot(point_t,point(:,p1_SWv,p),'g ')
plot(point_t,point(:,p1_Swnet,p),'k')
plot(point_t,point(:,p1_Swup,p),'r')
legend('Swin ab. veg.', 'SWv ads. veg.', 'Swnet soil', 'Swup ab. veg')
title('SW rad [W/m^2]','fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks');
%xlim([tbeg_pnt tend_pnt]);
% xlim([p_start p_end]);
% datetick('x','dd-mmm HH','keeplimits'); 
dynamicDateTicks
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 13: plot SWin partitioning'
%%%
% p1_Swin         =	22	%	[W/m2]
% p1_Swbeam       =	23	%	[W/m2]
% p1_Swdiff       =	24	%	[W/m2]
figure(13),clf
plot(point_t,point(:,p1_Swin,p)), hold on
plot(point_t,point(:,p1_Swbeam,p),'c ')
plot(point_t,point(:,p1_Swdiff,p),'r')
legend('SWin', 'SWin beam','SWin diff')
title('SWin partitioning [W/m^2]','fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks');
%xlim([tbeg_pnt tend_pnt]);
%datetick('x','dd/mm/yy','keeplimits'); %diplays the Months, Days and Years
% xlim([p_start p_end]);
% datetick('x','dd-mmm HH','keeplimits'); 
dynamicDateTicks

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp 'Fig 14: plot LW rad'
%%%
%  p1_Lwin p1_LWv p1_Lwnet p1_Lwup
figure(141),clf
plot(point_t,point(:,p1_Lwin,p)), hold on
plot(point_t,point(:,p1_LWv,p),'g ')
plot(point_t,point(:,p1_Lwnet,p),'k')
plot(point_t,point(:,p1_Lwup,p),'r')
plot(point_t,point(:,p1_Lwin,p)-point(:,p1_Lwup,p),'m')
legend('LWin', 'LWv ads.veg.', 'Lwnet soil', 'Lwup ab.veg.', 'Lwin - Lwup')
title('LW rad [W/m^2]','fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks');
%xlim([p_start p_end]);
%xlim([tbeg_pnt tend_pnt]);
%datetick('x','dd-mmm HH','keeplimits'); 
dynamicDateTicks

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig15: plot  energy balance at the soil/snow surface'
%%%
% p1_Swnet        =	28	%	[W/m2]
% p1_Lwnet        =	29	%	[W/m2]
% p1_H            =	30	%	[W/m2]
% p1_LE           =	31	%	[W/m2]
% p1_EB_surf      =	20	%	[W/m2]
Rn_net=point(:,p1_Swnet,p)+point(:,p1_Lwnet,p);

figure(15),clf
plot(point_t,Rn_net,'y'), hold on
plot(point_t,point(:,p1_H,p),'r')
plot(point_t,point(:,p1_LE,p),'b')
plot(point_t,point(:,p1_EB_surf,p),'g')
legend('Swnet+LWnet','H','LE','EB surf')%,'Qrain[W/m2]')
title('Energy balance at the soil/snow surface [W/m^2]','fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig15.1: plot  soil heat fluxes'
%%%
% p1_EB_surf      =	20	%	[W/m2]
% p1_Soil_heat_fl =	21	%	[W/m2]

figure(151),clf
plot(point_t,point(:,p1_EB_surf,p),'g'), hold on
plot(point_t,point(:,p1_Soil_heat_fl,p),'k')
legend('EB surf','Soil heat fl')%,'Qrain[W/m2]')
title('Soil heat fluxes [W/m^2]','fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 15.2: plot sensible fluxes components'
%%%
% p1_Hg_unveg     =	41	%	[W/m2]
% p1_Hg_veg       =	43	%	[W/m2]
% p1_Hv           =	39	%	[W/m2]
% p1_Estored_can	=	36	%	[W/m2]

figure(152),clf
plot(point_t,point(:,p1_Hg_unveg,p),'r'), hold on
plot(point_t,point(:,p1_Hg_veg,p),'m')
plot(point_t,point(:,p1_Hv,p),'g')
plot(point_t,point(:,p1_Estored_can,p),'k')
legend('Hg unveg', 'Hg veg', 'Hv' , 'Estored can')
title('sensible fluxes components [W/m^2]','fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 15.3: plot latent fluxes components'
%%%
% p1_LEg_unveg	=	42	%	[W/m2]
% p1_LEg_veg      =	44	%	[W/m2]
% p1_LEv          =	40	%	[W/m2]
% p1_Estored_can	=	36	%	[W/m2]

figure(153),clf
plot(point_t,point(:,p1_LEg_unveg,p),'r'), hold on
plot(point_t,point(:,p1_LEg_veg,p),'m')
plot(point_t,point(:,p1_LEv,p),'g')
plot(point_t,point(:,p1_Estored_can,p),'k')
legend('LEg unveg', 'LEg veg', 'LEv' , 'Estored can')
title('latent fluxes components [W/m^2]','fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 15.4: plot latent fluxes components parallel calculation'
%%%
%         Recalculate LE_up using the parallel resistance approach
%         LEv=(qv-qca)/rv   ->rv
%         LEgv=(qs-qca)/rg   -> rg
%         ra=(rv*rg)/(rv+rg) -> ra
%         LEup_veg=(qca-qa)/ra -> LEup_veg
%         LEup=fc*LEup_veg+(1-fc)*LEg

% LEup
% LEup_veg
% p1_LEg_unveg	=	42	%	[W/m2]
% p1_LEg_veg      =	44	%	[W/m2]
% p1_LEv          =	40	%	[W/m2]
% p1_Estored_can	=	36	%	[W/m2]

figure(154),clf
plot(point_t,LEup_veg(:,p),'g'), hold on
plot(point_t,LEup(:,p),'r')
% plot(point_t,point(:,p1_LEg_veg,p),'m')
% plot(point_t,point(:,p1_LEv,p),'g')
legend('LEup veg','LEup', 'LEg veg', 'LEv' , 'Estored can')
title('latent fluxes parallel calculation [W/m^2]','fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 16: plot total energy balance above canopy'
%%%
%Sw_net_atm=p1_Swin - p1_Swup
Sw_net_atm=point(:,p1_Swin,p) - point(:,p1_Swup,p);
%Lw_net_atm=p1_Lwin - p1_Lwup
Lw_net_atm=point(:,p1_Lwin,p) - point(:,p1_Lwup,p);
Rn_net_atm=Sw_net_atm+Lw_net_atm;
%p1_Hup    
%p1_LEup
%EB_ab_veg=Rn_net_atm-p1_Hup-p1_LEup
EB_ab_veg=Rn_net_atm-point(:,p1_Hup,p)-point(:,p1_LEup,p);

figure(16),clf
plot(point_t,Rn_net_atm,'y'), hold on
plot(point_t,EB_ab_veg,'g')
plot(point_t,point(:,p1_Hup,p),'r')
plot(point_t,point(:,p1_LEup,p),'b')
legend('Rnet tot', 'EB.ab.veg.', 'Hup' , 'LEup')
title('total energy balance above canopy [W/m^2]','fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 17: plot cumulated energy balance above canopy'
%%%
% cumulated energy balance 
EB_ab_veg_cum=cumsum(EB_ab_veg);
Rnet_atm_cum=cumsum(Rn_net_atm);
Hup_cum=cumsum(point(:,p1_Hup,p));
LEup_cum=cumsum(point(:,p1_LEup,p));

figure(17),clf
disp 'plot cumulated energy balance above canopy W m^{-2} h'
plot(point_t,Rnet_atm_cum,'m'), hold on
plot(point_t,Hup_cum,'r ')
plot(point_t,LEup_cum,'b')
plot(point_t,EB_ab_veg_cum,'g')
legend('Rnet.atm_{cm}','Hup_{cm}','LEup_{cm}','EB.ab.veg_{cm}'...
    ,'location','northwest')
title('cumulated energy balance above canopy [MJ]','fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 17.1: plot cumulated energy balance snow/soil surface'
%%%
% p1_Swnet        =	28	%	[W/m2]
% p1_Lwnet        =	29	%	[W/m2]
% p1_H            =	30	%	[W/m2]
% p1_LE           =	31	%	[W/m2]
% p1_EB_surf      =	20	%	[W/m2]
%Rn_net=point(:,p1_Swnet)+point(:,p1_Lwnet);

% cumulated energy balance 
EB_surf_cum=cumsum(point(:,p1_EB_surf,p));
Rn_net_cum=cumsum(Rn_net);
H_cum=cumsum(point(:,p1_H,p));
LE_cum=cumsum(point(:,p1_LE,p));

figure(171),clf
disp 'plot cumulated energy balance snow/soil surface W m^{-2} h'
plot(point_t,Rn_net_cum,'m'), hold on
plot(point_t,H_cum,'r ')
plot(point_t,LE_cum,'b')
plot(point_t,EB_surf_cum,'g')
legend('Rnet_{cm}','H_{cm}','LE_{cm}','EB.surf_{cm}'...
    ,'location','northwest')
title('cumulated energy balance snow/soil surface [MJ]','fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 18: plot  evaporation partitioning'
%%%
% p1_Evap_surf	=	45	%	mm
% p1_Trasp_can	=	46	%	mm' 
figure(18),clf
plot(point_t,point(:,p1_Evap_surf,p),'b'), hold on
plot(point_t,point(:,p1_Trasp_can,p),'g')
plot(point_t,(point(:,p1_Trasp_can,p))+(point(:,p1_Evap_surf,p)),'k')
legend('Evap surf','Trasp can','ET total')
title([' evaporation partitioning [mm/Dt]. Totals Eg = ',...
    num2str(floor(sum(point(:,p1_Evap_surf,p)))),' mm, Etc = ',...
    num2str(floor(sum(point(:,p1_Trasp_can,p)))),' mm, Et tot = ',...
    num2str(floor(sum(point(:,p1_Trasp_can,p))+sum(point(:,p1_Evap_surf,p)))),' mm'...
    ],'fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 181: plot  cumulated evaporation partitioning'
%%%



figure(181),clf
plot(point_t,cumsum(point(:,p1_Evap_surf,p)),'b'), hold on
plot(point_t,cumsum(point(:,p1_Trasp_can,p)),'g')
plot(point_t,cumsum(point(:,p1_Trasp_can,p))+cumsum(point(:,p1_Evap_surf,p)),'k')
%plot(point_t,point(:,46),'g')
legend('Evap surf','Trasp can','ET total','location','northwest')
title([' cumulated evaporation partitioning [mm]. Totals Eg = ',...
    num2str(floor(sum(point(:,p1_Evap_surf,p)))),' mm, Etc = ',...
    num2str(floor(sum(point(:,p1_Trasp_can,p)))),' mm'],'fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks

%% cumulated check - energy and mass budget
LEup_mm=ET_conversion(point(:,p1_LEup,p),point(:,p1_Tair,p),3600);
LE_mm=ET_conversion(point(:,p1_LE,p),point(:,p1_Tair,p),3600);
LEv_mm=ET_conversion(point(:,p1_LEv,p),point(:,p1_Tair,p),3600);
LEg_veg_mm=ET_conversion(point(:,p1_LEg_veg,p),point(:,p1_Tair,p),3600);
LEg_unveg_mm=ET_conversion(point(:,p1_LEg_unveg,p),point(:,p1_Tair,p),3600);

LEup_mm_c=cumsum(LEup_mm);
LE_mm_c=cumsum(LE_mm);
LEv_mm_c=cumsum(LEv_mm);
LEg_veg_mm_c=cumsum(LEg_veg_mm);
LEg_unveg_mm_c=cumsum(LEg_unveg_mm);

% plot check - energy and mass budget
figure(182),clf
plot(point_t,cumsum(point(:,p1_Evap_surf,p)),'r','LineWidth',2), hold on
plot(point_t,cumsum(point(:,p1_Trasp_can,p)),'g','LineWidth',2)
plot(point_t,cumsum(point(:,p1_Trasp_can,p))+cumsum(point(:,p1_Evap_surf,p)),'k','LineWidth',2)
plot(point_t,LEup_mm_c,'--k','LineWidth',2)
plot(point_t,LE_mm_c,'k')
plot(point_t,LEv_mm_c.*point(:,p1_fc,p),'--g','LineWidth',2)
plot(point_t,LEg_veg_mm_c.*point(:,p1_fc,p),'--b','LineWidth',2)
plot(point_t,LEg_unveg_mm_c.*(1-point(:,p1_fc,p)),'m','LineWidth',2)
%plot(point_t,point(:,46),'g')
legend('Evap surf','Trasp can','ET total',...
    'LEup','LE','LEv*fc','LEg veg * fc','LEg unveg *(1-fc)','location','northwest')
title([' cumulated evaporation partitioning [mm]. Totals Eg = ',...
    num2str(floor(sum(point(:,p1_Evap_surf,p)))),' mm, Etc = ',...
    num2str(floor(sum(point(:,p1_Trasp_can,p)))),' mm'],'fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks

%% cumulated check - energy and mass budget
% alternative approach using parallel resitance approach
%         LEv=(qv-qca)/rv   ->rv
%         LEgv=(qs-qca)/rg   -> rg
%         ra=(rv*rg)/(rv+rg) -> ra
%         LEup_veg=(qca-qa)/ra -> LEup_veg
%         LEup=fc*LEup_veg+(1-fc)*LEg

LEup1_mm=ET_conversion(LEup(:,p),point(:,p1_Tair,p),3600);

LEup1_mm_c=cumsum(LEup1_mm);

% plot check - energy and mass budget
figure(183),clf
plot(point_t,cumsum(point(:,p1_Evap_surf,p)),'r','LineWidth',2), hold on
plot(point_t,cumsum(point(:,p1_Trasp_can,p)),'g','LineWidth',2)
plot(point_t,cumsum(point(:,p1_Trasp_can,p))+cumsum(point(:,p1_Evap_surf,p)),'k','LineWidth',2)
plot(point_t,LEup1_mm_c,'--','color',[0.5 0.5 0.5],'LineWidth',2)
plot(point_t,LE_mm_c,'k')
plot(point_t,LEv_mm_c.*point(:,p1_fc,p),'--g','LineWidth',2)
plot(point_t,LEg_veg_mm_c.*point(:,p1_fc,p),'--b','LineWidth',2)
plot(point_t,LEg_unveg_mm_c.*(1-point(:,p1_fc,p)),'m','LineWidth',2)
%plot(point_t,point(:,46),'g')
legend('Evap surf','Trasp can','ET total',...
    'LEup1 parallel','LE','LEv*fc','LEg veg * fc','LEg unveg *(1-fc)','location','northwest')
title([' cumulated ET partitioning - parallel approach [mm]. Totals Eg = ',...
    num2str(floor(sum(point(:,p1_Evap_surf,p)))),' mm, Etc = ',...
    num2str(floor(sum(point(:,p1_Trasp_can,p)))),' mm'],'fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 19: plot  precipitation interception '
%%%

figure(19),clf
plot(point_t,point(:,p1_Water_on_can,p),'r'), hold on
plot(point_t,point(:,p1_Snow_on_can,p),'b')
legend('Wtrain[mm] ','Wtsnow[mm]')
title([' precipitation interception [mm/Dt]'],'fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 20: plot below canopy precipitation partitioning '
%%%

figure(20),clf
plot(point_t,point(:,p1_Prain_bel_can,p),'b'), hold on
plot(point_t,point(:,p1_Psnow_bel_can,p),'r')
plot(point_t,point(:,p1_Prain_on_snow,p),'--k')
legend('P_{rain}','P_{snow}','P_{Prain on snow}')
title([' below canopy precipitation partitioning [mm/Dt]'],'fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 21: plot below canopy  cumulated precipitation components '
%%%
Prain_bel_can_cum=cumsum(point(:,p1_Prain_bel_can,p));
Psnow_bel_can_cum=cumsum(point(:,p1_Psnow_bel_can,p));
Ptot_bel_can_cum=Prain_bel_can_cum+Psnow_bel_can_cum;

figure(21),clf
plot(point_t,Prain_bel_can_cum,'b'), hold on
plot(point_t,Psnow_bel_can_cum,'r')
plot(point_t,Ptot_bel_can_cum,'k')
legend('Prain c[mm]','Psnow c[mm]','Ptot c[mm]',...
    'location','northwest')
title(' below canopy cumulated precipitation components [mm]','fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 22: plot above canopy precipitation partitioning '
%%%

figure(22),clf
plot(point_t,point(:,p1_Prain_abo_can,p),'b'), hold on
plot(point_t,point(:,p1_Psnow_abo_can,p),'r')
legend('P_{rain}','P_{snow}')
title([' above canopy precipitation partitioning [mm/Dt]'],'fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 23: plot above canopy cumulated precipitation partitioning '
%%%
Prain_abo_can_cum=cumsum(point(:,p1_Prain_abo_can,p));
Psnow_abo_can_cum=cumsum(point(:,p1_Psnow_abo_can,p));
Ptot_abo_can_cum=Prain_abo_can_cum+Psnow_abo_can_cum;

figure(23),clf
plot(point_t,Prain_abo_can_cum,'b'), hold on
plot(point_t,Psnow_abo_can_cum,'r')
plot(point_t,Ptot_abo_can_cum,'k')
legend('Prain c[mm]','Psnow c[mm]','Ptot c[mm]',...
    'location','northwest')
title([' above canopy  cumulated precipitation partitioning [mm]'],'fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 24: plot snow melt '
%%%
%     p1_Snow_melt	=	65;	%	mm
%     p1_Snow_subl	=	66;	%	mm
%     p1_Snow_blown	=	67;	%	mm
figure(24),clf
plot(point_t,point(:,p1_Snow_melt,p),'b'), hold on
plot(point_t,point(:,p1_Snow_subl,p),'r')
plot(point_t,point(:,p1_Snow_blown,p),'k')
legend('Snow melt','Snow evap','Snow subl')
title([' snow melt [mm/Dt]'],'fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 25: plot cumulated snow melt '
%%%

Snow_melt_cum=cumsum(point(:,p1_Snow_melt,p));
Snow_subl_cum=cumsum(point(:,p1_Snow_subl,p));
Snow_blown_cum=cumsum(point(:,p1_Snow_blown,p));
figure(25),clf
plot(point_t,Snow_melt_cum,'b'), hold on
plot(point_t,Snow_subl_cum,'r')
plot(point_t,Snow_blown_cum,'k')
legend('Snow melt c','Snow evap c','Snow subl c')
title(['Cumulated snow melt [mm]'],'fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 26: plot water tables '
%%%
% p1_depth_thaw	=	65	%	mm
% p1_depth_wat	=	66	%	mm
figure(26),clf
plot(point_t,-point(:,p1_depth_thaw,p),'b'), hold on
plot(point_t,-point(:,p1_depth_wat,p),'k')
legend('thawed soil depth','water table depth')
title([' water tables [mm]'],'fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp 'Fig 27: plot SWE and snow depth '
%%%
%     p1_Snow_depth = 61;  %   mm
%     p1_Snow_swe = 62;  %   mm
figure(27),clf
% snow depth
subplot(2,1,1) 
plot(point_t,point(:,p1_Snow_depth,p)/10,'b')
legend('Snow depth [cm]')
title([' Snow depth [cm] '],'fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks
% SWE
subplot(2,1,2) 
plot(point_t,point(:,p1_Snow_swe,p),'b')
legend('SWE [mm]')
title([' SWE [mm] '],'fontsize',12)
set(gca,'fontsize',12), grid on
datetick('x',12,'keepticks'), %xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp 'Fig 20: plot  runoff '
% %%%
% %  62)Runoff[mm] 63)'q_sup[mm]', 64)'q_sub[mm]' 65)'DS_sup[mm]',
% %  66)'DS_sub[mm]', 67)'q_G[mm]'
% figure(20),clf
% plot(point_t,point(:,62),'b'), hold on
% plot(point_t,point(:,63),'g')
% plot(point_t,point(:,64),'m')
% plot(point_t,point(:,65),'--b')
% plot(point_t,point(:,66),'--m')
% plot(point_t,point(:,67),'k')
% legend('62 Runoff [mm/Dt]','63 q sup[mm]', '64 q sub[mm]','65 DS sup[mm]','66 DS sub[mm]','67 q G[mm]')
% title('Runoff components [mm/Dt]','fontsize',12)
% set(gca,'fontsize',12), grid on
% datetick('x',12,'keepticks');

