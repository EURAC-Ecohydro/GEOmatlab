%% plot_GEOtop_2_0_profiles.m
% Copyright, 2010 Giacomo Bertoldi
%
% plot 2_0 GEOTOP point profiles for every layer
% of: soil temparature, sil pressure, soil total water content, soil ice
% content
%
% load before load_GEOtop_2_0
%%

% define here if you have multiple point output the point to plot
mypoint=1


mylegend=num2str([1:Nz]');
disp(sprintf('Simulation %s',mysim));
disp(sprintf('Ouptut plotted on %s',date));

if(ifbasin==1)
    tbeg_bas=min(basin_t); tend_bas=max(basin_t);
end
tbeg_pnt=min(point_t); tend_pnt=max(point_t);

%%
% total soil water [mm]
vliq_tot=thetaz;
for n=1:npoints
    for k=1:Nz
        vliq_tot(:,noffset+k,n)=(thetaz(:,noffset+k,n))*Dz(k);
    end
end

% cumulative total soil water [mm]
vliq_cum=vliq_tot;
tmp=vliq_tot(:,noffset+1:noffset+Nz,:);
tmp_cum=cumsum(tmp,2);
vliq_cum(:,noffset+1:noffset+Nz,:)=tmp_cum;
clear tmp tmp_cum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOIL TEMPERATURE PROFILE
% Date12[DDMMY6YYhhmm], 1) JulianDayFromYear0[days], 2) TimeFromStart[days], 3) Simulation_Period,Run, 4) IDpoint,
% 5)'Layer1 ' 4)'Layer2 ' 7)'... 'Layer Nz'
disp 'Tz file : T Layer 1..Nz in Tz[C]'
figure(30),clf
disp 'plot soil temperatues [C]'
plot(Tz_t,Tz(:,noffset+1:Nz+noffset,mypoint)), hold on
% 95)'Tg[C]=ground or snow surface '
plot(point_t,point(:,p1_Tair,mypoint),'--k')
title('soil temperatues [C]','fontsize',12)
xlabel('t [days]','fontsize',12)
set(gca,'fontsize',12), grid on
legend(mylegend)
datetick('x',12,'keepticks'), xlim([tbeg_pnt tend_pnt]);
dynamicDateTicks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOIL MOISTURE PROFILE
% Date12[DDMMY6YYhhmm], 1) JulianDayFromYear0[days], 2) TimeFromStart[days], 3) Simulation_Period,Run, 4) IDpoint,
% 5)'Layer1 ' 4)'Layer2 ' 7)'... 'Layer Nz'
if (ifthetaz)
    disp 'po06_THETAz0001.txt in thetaz'
    figure(31),clf
    disp 'plot soil \theta [ ]'
    plot(thetaz_t,thetaz(:,noffset+1:Nz+noffset,mypoint))
    title('soil \theta [ ]','fontsize',12)
    xlabel('t [days]','fontsize',12)
    set(gca,'fontsize',12), grid on
    legend(mylegend)
    datetick('x',12,'keepticks')% xlim([tbeg_pnt tend_pnt]);
    dynamicDateTicks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% RESCALED SOIL MOISTURE PROFILE
    % Date12[DDMMY6YYhhmm], 1) JulianDayFromYear0[days], 2) TimeFromStart[days], 3) Simulation_Period,Run, 4) IDpoint,
    % 5)'Layer1 ' 4)'Layer2 ' 7)'... 'Layer Nz'
    
    
    disp 'po06_THETAz0001.txt in thetaz'
    figure(311),clf
    disp 'plot rescaled soil \theta [ ]'
    %plot(thetaz_t,satz(:,3)), hold on
    plot(thetaz_t,satz(:,noffset+1:Nz+noffset,mypoint))
    title('soil saturation=(\theta-\theta_r)/(\theta_s-\theta_r) [ ]','fontsize',12)
    xlabel('t [days]','fontsize',12)
    set(gca,'fontsize',12), grid on
    legend(mylegend)
    datetick('x',12,'keepticks')% xlim([tbeg_pnt tend_pnt]);
    dynamicDateTicks
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WATER STORAGE PROFILE
    
    disp 'po06_THETAz0001.txt in thetaz'
    figure(312),clf
    disp 'plot soil \theta [ ]'
    plot(thetaz_t,vliq_tot(:,noffset+1:Nz+noffset,mypoint))
    title('soil liquid water content each layer [mm]','fontsize',12)
    xlabel('t [days]','fontsize',12)
    set(gca,'fontsize',12), grid on
    legend(mylegend)
    datetick('x',12,'keepticks')% xlim([tbeg_pnt tend_pnt]);
    dynamicDateTicks
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WATER STORAGE CUMULATED PROFILE
        
    disp 'po06_THETAz0001.txt in thetaz'
    figure(313),clf
    disp 'plot soil \theta [ ]'
    plot(thetaz_t,vliq_cum(:,noffset+1:6+noffset,mypoint))
    title('cumulative soil liquid water content each layer [mm]','fontsize',12)
    xlabel('t [days]','fontsize',12)
    set(gca,'fontsize',12), grid on
    legend(mylegend)
    datetick('x',12,'keepticks')% xlim([tbeg_pnt tend_pnt]);
    dynamicDateTicks
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOIL PRESSURE PROFILE
% Date12[DDMMY6YYhhmm], 1) JulianDayFromYear0[days], 2) TimeFromStart[days], 3) Simulation_Period,Run, 4) IDpoint,
% 5)'Layer1 ' 4)'Layer2 ' 7)'... 'Layer Nz'

% correct psi from mm in m and eliminate positive values to plot in log scale
if (ifpsiz)
    mypsiz=psiz/1000;
    tmp=find(mypsiz>=0);
    mypsiz(tmp)=NaN;
    
    disp 'po05_psiz0001.txt.txt in psiz'
    figure(32),clf
    disp 'plot soil \psi [m]'
    semilogy(psiz_t,mypsiz(:,noffset+1:Nz+noffset,mypoint))
    title('soil soil \psi [m]','fontsize',12)
    xlabel('t [days]','fontsize',12)
    set(gca,'fontsize',14), grid on
    legend(mylegend)
    datetick('x',12,'keepticks'), xlim([tbeg_pnt tend_pnt]);
    %ylim([-200 1]);
    dynamicDateTicks
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOIL ICE CONTENT PROFILE
% Date12[DDMMY6YYhhmm], 1) JulianDayFromYear0[days], 2) TimeFromStart[days], 3) Simulation_Period,Run, 4) IDpoint,
% 5)'Layer1 ' 4)'Layer2 ' 7)'... 'Layer Nz'
if (ifthetaicez)
    disp '08_THETAICEz0001.txt in thetaicez'
    figure(33),clf
    disp 'plot ice \theta [ ]'
    plot(thetaz_t,thetaicez(:,noffset+1:Nz+noffset,mypoint))
    title('ice \theta [ ]','fontsize',12)
    xlabel('t [days]','fontsize',12)
    set(gca,'fontsize',12), grid on
    legend(mylegend)
    datetick('x',12,'keepticks'), xlim([tbeg_pnt tend_pnt]);
    dynamicDateTicks
end