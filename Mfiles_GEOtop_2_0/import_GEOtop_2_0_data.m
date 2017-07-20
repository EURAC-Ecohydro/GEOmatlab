function [mydata, mytimes, mydates, mylabels]=import_GEOtop_2_0_data(myfile)

%function [mydata, mytimes, mydates, mylabels]=import_GEOtop_2_0_data(myfile)
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

% column with time data
%coltime=1; % GEOtop 1_32
%coltime=2; % GEOtop 1_1
coltime=1; % GEOtop 1_2 and 2_0

tmp=importdata(myfile,',',1);

    s=whos('tmp');   
    if (strcmp(s.class,'struct'))
        mydata=tmp.data;
        mydates=tmp.textdata(2:size(tmp.textdata,1),1);
        mytimes=datenum(mydata(:,coltime));
        %mytimes=datenum(mydates,'dd/mm/yyyy HH:MM')/1000
        mylabels=tmp.textdata(1,2:size(tmp.textdata,2));
        tbeg=min(mytimes);
        clear tmp
    else
        disp(fprintf('\nWarning: no data in %s \n',myfile))
        mydata=0, mytimes=0, mydates=0, mylabels='';
        clear tmp
    end
    
   
return
