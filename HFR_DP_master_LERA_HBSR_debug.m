%HFR_DP_master_XXXX.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HFR_DP_master_XXXX.m
%
%  This script contains the master program that runs the HFR_DP program of
%  spectra to radial average processing for LERA-type HF radar datasets.
%
%  SETUP:  
%   To setup this template file or a particular radar site. 
%  The user must:
%  (1) define the site, array, directory locations, and processing 
%      choices in the first portion of the script, 
%  (2) save a new version of this template file with the XXXX replaced by
%      the site name
%  (3) Run the script within matlab, from the command line, or via a cron job to process
%      timeseries files to radial averaged files.
%
%  See HFR_DP_SETUP_README.m for directions on how to set up a machine to
%  perform the processing
%
% v1        March 2017
%
% Anthony Kirincich
% Woods Hole Oceanographic Institution
% akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%

clear
%%% load hamlin paths for HFR_Progs, m_map, etc.
%cd /home/hfrp/ ; HFR_DP_addpaths_hamlin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     begin user choices       %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set site/array information

%%% the site name to process
site_name='HBSR';
Site_loc = [41+30.5385/60 -(70+4.1745/60)];
Site_bounds=[90 90+180];

%site_name='NWTP';
%Site_loc = [41+14/60+30.96/3600 -(70+6/60+24.85/3600) ];
%Site_bounds=[120 120+180];
    
%%% set the common constants for the HFR array this site is a part of:
ARRAY=[];
ARRAY.name='NES';  
%%% array name, does not have to be four letters like radial sites, can be the
%%% same name for multiple radar sites that contribute to coverage within an area

%%% define the bounding box of the array, needed for any plotting operations.
ARRAY.min_lon=-(71+50/60); 
ARRAY.max_lon=-(69+40/60); 
ARRAY.min_lat=(40+40/60); 
ARRAY.max_lat=(41+40/60);

%%% Note that if you change the map data or the bounding box of the ARRAY
%%% without changing the ARRAY name, you will need to remove/delete the
%%% 'ARRAY.name'_coast.mat file from the working directory.  
%%%
%%% The mapping script only (for speed) regenerates the coastline if there
%%% is no coastline file ('ARRAY.name'_coast.mat) found.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set directory choices

%%% What is the working directory where the HFR_DP files are located
working_dir='/Users/anthony/Matlab/working/lera/lera_DPd';

%%% Where are the scripts
scripts_dir=[working_dir '/scripts/'];

%%% Where are the data  (see HFR_DP_SETUP_README.m for instructions
base_dir='/Users/anthony/Data/LERA_process';
%base_dir='/Codar/SeaSonde/Data/RadialSites'; %if wish to place in codar directory
%base_dir='/Volumes/hfrp/bdata/RadialSites';
%base_dir='/Volumes/data/RadialSites';

%%% Where do the DP logs go
log_dir=['/Users/anthony/Matlab/working/LOGS'];
%log_dir='/Codar/SeaSonde/Logs/codar_DP_Logs'; %if wish to place in codar directory


eval(['addpath ' scripts_dir]);
eval(['addpath ' working_dir]);
%%% move to inital base directory
eval(['cd ' working_dir])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set constants for the HFR
CONST=[];
%%% constants that won't change %%%
CONST.c=299792458;  %m/s speed of light in a vaccum
    %    c=299705000;   % in air   
    %    c=299644310;   % in air based on getting v=0 at ibragg (difference of 2 mm/s from first value....)
%number of antenna elements
CONST.N=8;
%%% how many solutions to process
CONST.Nmax=4;

%%% set info on the antenna makeup
CONST.RxAntConfig='8-channel Rectangular Array';
CONST.TxAntConfig='4-post Quad Array';
CONST.Tx_bearing= '180';  %in degT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set parameters for spectral processing

%%% the final azimuthal resolution of the radial velocities
CONST.bearing_width=5;

%%% the minimum Signal to Noise ratio allowed (used in Radial averaging and ImageFOLs     
CONST.snr_min=5;  %dB

%%%% the endpoints of the ave or RM pattern used for DF
% % NWTP
%CONST.pattern_ends=[90 276];  % in degT
%i=find(patt.angles>90 & patt.angles<276);
%LPWR
CONST.pattern_ends=[90 270];  % in degT


%%% for imageFOL, set user parameters
%%%  parameters are: [velD_change max_vel snr_min];  with vels defined in cm/s
%%%  see imageFOL for more details
%CONST.imageFOL_user_param=[15 100 CONST.snr_min];  %for a 25MHz site
%CONST.imageFOL_user_param=[20 100 CONST.snr_min];  %for a 25MHz site
%CONST.imageFOL_user_param=[40 300 CONST.snr_min];  %for a 5MHz site, viewing the Gulf Stream
CONST.imageFOL_user_param=[20 100 CONST.snr_min];  %for a 16MHz site
%CONST.imageFOL_user_param=[40 100 CONST.snr_min];  %for a 16MHz site

%%% define radave RM thresholds... to be used to weed out bad radials before the radial average is made 
%CONST.radave_thresholds=[snr_thresh angpeak_thresh angwidth_thresh(1) angwidth_thresh(2)];
%CONST.radave_thresholds=[CONST.snr_min 5 0 50];
CONST.radave_thresholds=[CONST.snr_min .001 0 200]; % for mle or wsf


%%% define the peak threshold to use for the music doa peak finder,
%%% this is done in 10*log10(MS) space, making .05-.5 reasonable thresholds
%%% that scale with the values to allow small peaks to exist.
%CONST.doa_peak_thresh=.05;
CONST.doa_peak_thresh=.25;
%CONST.doa_peak_thresh=.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set processing choices

%%% Which first order line method to use
CONST.FOL_type='imageFOL';    %uses Image FOL methods
%CONST.FOL_type='codarFOL';    %imports the COS FOLs from the spectral file

%%% Which radial averaging method to output, both are estimated in the
%%% radave step, but only one can be sent out within the lluv file format.
CONST.radave_type='regular';    %follows COS arthmetic averaging of the radials
                                %  within the CONST.bearing_width area and error calcuation
%CONST.radave_type='metricQC';   %follows Kirincich et al 2012 to use thresholding
                                %  and power-weighted spatial means for rad ave

%%% which type of antenna manifold will be used
CONST.which_patt='ideal';                                
%CONST.which_patt='meas';                                
              
%%% which type of direction finding will be used
CONST.which_df='music';                                
%CONST.which_df='wsf';                                
%CONST.which_df='mle';                                

%%%% which type of method to determine the number of emitters...see Johnson and Dudgen, sec 7.3.5 for details
%CONST.which_Ns_meth='MDL';
%CONST.which_Ns_meth='AIC';
%CONST.which_Ns_meth='music_param';
CONST.which_Ns_meth='music_highest';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set processing constants, the only user selectable parts
CONST.PC=[];
CONST.PC.MaxRangeKM=100;        % the maximum ranges, in distance to be kept in spectra, from >0 to Max
%CONST.PC.SpecOverlapPct=50;     %for WOSA, style increased ensembles into CSE

%CONST.PC.SpeclengthPnts=2048;   %defines the number of points in the spectra
CONST.PC.SpecOverlapPct=78;     %for WOSA, style increased ensembles into CSE
CONST.PC.SpeclengthPnts=1024*2;   %defines the number of points in the spectra
%CONST.PC.SpeclengthPnts=512;   %defines the number of points in the spectra



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set plotting and printing switches
%plotting:
CONST.goplot=[1 1];   %where:
%  switch 1 controls the final plots of radial results, etc.  (1=plot, 0=do not plot)
%  switch 2 controls the intermediary plots within functions  (1=plot, 0=do not plot)
CONST.goprint=[0 0];   
%  switch 1 prints the final plots of radial results, etc.  (1=print, 0=do not print)
%  switch 2 prints the intermediary plots within functions  (1=print, 0=do not print)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set how 'new' files to process are found

%%% Method 1: choose files based on the file time
CONST.files_to_process_method=1;
%%% Method 2: compare css files to RM files, and process new files
%CONST.files_to_process_method=2;

%%% set the time frame to examine, only used by method 1
% start_time=datenum(2017,8,15,12,0,0);
% %start_time=datenum(2017,8,17,0,00,0);
% end_time=datenum(2017,8,18,0,0,0);

start_time=datenum(2018,11,2,0,0,0);
end_time=datenum(2018,12,1,0,0,0);
CONST.files_to_process_dates=[start_time end_time];


%%%%%%%%%%% end radial processing setup %%%%%%%%%%    
%%    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% running the processing scripts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%% setup logfile name for diary
%  logfn=['LOG_v7_' mfilename '_'  datestr(now,30) '.log'];
%  diary([log_dir '/' logfn]);
%  tic

%%% for this site, prepare for processing
%HFR_DP_lera_ts2radial_prepwork_v2  
HFR_DP_lera_ts2radial_prepwork_v3
%%% sets up processing directories, loads site info and cal, and 
%%% identifies files that need processing

%  %limit the number of files  on in any 1 pass of this processing.  
%  % this can be helpful if running in realtime
%  if length(fnames)>5;
%      fnames=fnames(1:5);
%  end

%%
%%% now loop overeach file in fnames to process data
for jjj=1:1:length(fnames)
    
%%% display the file to be processed
jjj=13
filein=fnames{jjj}

%%
disp(filein)

% %reload header and patt files
% eval(['cd ' incoming_config_file_dir])
% %load site header
% config_dir=[incoming_config_file_dir '/'];
% radar_header %radar_header is used here as well as on the site computer to load this info for the compression step
% % thus don't change it now.
% eval(['cd ' working_dir])

%%%% start processing steps list now that HEAD has been established.
HEAD.ProcessingSteps={mfilename};

%CONST.goplot(2)=1;
%%% load the ts data for this file and make the ensemble-averaged spectra
%[SpecHead,data,mtime]=lera_time2spectra_v2(incoming_ts_file_dir,incoming_spectra_file_dir,filein,RC,HEAD,CONST);
%[SpecHead,data,mtime]=lera_time2spectra_v4(incoming_ts_file_dir,incoming_spectra_file_dir,filein,RC,HEAD,CONST);
[SpecHead,data,mtime]=lera_time2spectra_v5(incoming_ts_file_dir,incoming_spectra_file_dir,filein,RC,HEAD,CONST);


%%%save the spectra file    
%%%% get the output filename, based on the filein %%%%
i=find(filein=='_');
%fileout=['CSE_' RC.SiteName '_' filein(i(2)+1:end)];
fileout=['CSE_' RC.SiteName '_' filein(1:i(1)-1) '.mat'];

%%%%% save this data for direction finding.
eval(['cd ' incoming_spectra_file_dir])
% save(fileout,'RC','SpecHead','data','mtime')
eval(['cd ' working_dir])
%%%%%%%%%%

%%% make a picture of the spectral estimate.
if CONST.goplot(1)==1;
    %set up plot of power for all ants
    figure(1); clf;
    r=SpecHead.RangeKm'*ones(1,length(SpecHead.doppler_velc));
    v=ones(length(SpecHead.RangeKm),1)*SpecHead.doppler_velc;
    
    for ii=1:RC.NANT
        eval(['d=data.a' num2str(ii) num2str(ii) ';'])
        subplot(4,2,ii);
        pcolor(v,r,20*log10(d)); shading flat; hg; caxis([-100 20]);
        plot([1; 1]*SpecHead.doppler_velc(SpecHead.iFBraggc),[1 1; max(r(:)) max(r(:))],'k')
        title(['Antenna ' num2str(ii) ' Ens. Averaged Autospectra'])
        if ii==7;  ylabel('range (km)'); xlabel('Doppler vel (m/s)'); end
        if ii==1; t=text(6,120,fileout); set(t,'interpreter','none'); end
    end
    c=colorbar;  set(c,'position',[.925 .1 .01 .2])
    
%     %%%% code to plot the ave power level at 2 ranges.
%     dd=[]; dd2=[];
%     figure(10); clf;
%      for ii=1:RC.NANT
%         eval(['d=data.a' num2str(ii) num2str(ii) ';'])
%         %subplot(4,2,ii);
%         dd(ii,:)=20*log10(d(10,:));
%         dd2(ii,:)=20*log10(d(20,:));
%    end
%    plot(v,mean(dd));  hg;
%    plot(v,mean(dd2));  hg;
%       axis([-8 8 -120 20])
%   title(fileout,'interpreter','none')
 
    
    if CONST.goprint(1)==1
        %print to figure
        figure(1); hold on
        set(gcf,'paperposition',[.25 .15 10 10])
        %%%%% save this data for direction finding.
        eval(['cd ' outgoing_pics_file_dir])
        i=find(fileout=='.');
        print('-djpeg','-r300',[fileout(1:i(1)-1) '_' date '.jpg'])
        eval(['cd ' working_dir])
    end
end

%%
%%%% if just wish to loop over spectral creation step
pause
end
return

%%
%%%% run radial metric processing %%%
%HFR_DP_lera_spectra2radialmetric_process_v1 

% for i5=1:3
% %%% which type of direction finding will be used
% if i5==1; CONST.which_df='music';                                
% elseif i5==2; CONST.which_df='wsf';                                
% elseif i5==3; CONST.which_df='mle';                                
% end

%HFR_DP_lera_spectra2radialmetric_process_v2   % use this one for DP_methods
%HFR_DP_lera_spectra2radialmetric_process_v3    %use this one for music v1
HFR_DP_lera_spectra2radialmetric_process_v4    %add data cut for near ranges, due to wave reflections

%%% if wishing to compare different RM for the same file
% eval(['RM' num2str(i5) '=RM;'])
% 

% %%%% if just doing the FOLS
% pause
%  end
 

%
%again, output format of RM(jjj).data is similar to COS RSv7 radial metric files
%with changes to the field contents to handle N=8
%mostly, only the results for one peak of one solution are given in
%each line,
% cols   fields
% 1-2    lon lat
% 3-4    u v   (here nan as will not be used)
% 5      flag    (here nan, used by COS but not here)
% 6 range
% 7 bearing
% 8 vel
% 9 direction
% 10  rangecell
% 11  dopcell
% 12  angselect (which solution if isave out is given
% 13  which peak of this solution is given
% 14 musicpow(v)
% 15 music_pkwidth
% 16 musicDOASpeak
% 17 spectraA3(snr)
% 18-20 musicEigenRatio musicpowerRatio musicoffRatio
% 21-27 which solutions were viable (1 yes 0 no, for 1-7 (need better description/name)
% 28-35 Eigenval (1-8)
% 36   DF_flag (see MUSIC algo above for explaination)        


%%% save this RM ouput...you don't have choice %%%
%move to saving directory
eval(['cd ' outgoing_radialmat_file_dir])
%f=char(fnames(jjj));
d=datevec(fdates(jjj));
%i=find(f=='_'); i2=find(f=='.');
%s=['RM_' site_name '_' num2str(d(1)) f(i(3):i2(1)) '.mat'];
MM=num2str(d(2)); if d(2)<10; MM=['0' MM]; end
dd=num2str(d(3)); if d(3)<10; dd=['0' dd]; end
hh=num2str(d(4)); if d(4)<10; hh=['0' hh]; end
mm=num2str(d(5)); if d(5)<10; mm=['0' mm]; end
file_date_str=[num2str(d(1)) '_' MM  '_' dd '_' hh mm];

s=['RM_' site_name '_' file_date_str '.mat'];
eval(['save ' s ' RM;'])
%move to inital base directory
eval(['cd ' base_dir])  

%
%%% plot the RM result on a map of the array domain %%%
if CONST.goplot(1)==1
    figure(3); clf;
    HFR_DP_quickmap(ARRAY,working_dir);
    title(['RM for Site:' RM.Site_name ', Time:' datestr(RM.time)])
    %m_plot(RM.data(:,1),RM.data(:,2),'k.')
    m_quiver(RM.data(:,1),RM.data(:,2),RM.data(:,3),RM.data(:,4),10,'r');
    
    %%% if saving figures.
    if CONST.goprint(1)==1;
        eval(['cd ' outgoing_pics_file_dir])
        
        %save pictures of spectra
        figure(7); hold on;
        i=find(s=='.');   s2=s(1:i(1)-1);
        t=text(-100,5,[s2 '      raw spectra w/ridge lines']); set(t,'interpreter','none','rotation',90);
        set(gcf,'paperposition',[.25 .25 10 6])
        print('-djpeg','-r300',[s2 '_sfols_' date '.jpg'])
        
        %save pictures of radials
        figure(3); hold on;
        t=title([s '      RM locations']); set(t,'interpreter','none');
        set(gcf,'paperposition',[.25 .25 5 5])
        print('-djpeg','-r300',[s2 '_rmlocs_' date '.jpg'])
        %move back to inital base directory
        eval(['cd ' working_dir])
    end% if go print
end% if goplot
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute Radial Average (equivalent to COS radial short result)

%%%% if RM is big enough continue to use RM to process radial metrics to radave files.
if length(RM.data)>500;
    
%%% compute radial averages   
[RADIAL,HEAD]=HFR_DP_lera_RadialAverage_v1(RM,patt,HEAD,CONST);

%%% Clean large radials using ImageFOL user parameter for max velocity found
for k = 1:numel(RADIAL)
    [Rclean,I]=cleanRadials_v7(RADIAL(k),CONST.imageFOL_user_param(2));
end

%%% save this Rclean ouput...you don't have choice %%%
%move to saving directory
eval(['cd ' outgoing_radialavemat_file_dir])
% f=char(fnames(jjj));
% d=datevec(fdates(jjj));
% i=find(f=='_'); i2=find(f=='.');
% s=['RadAve_' site_name '_' num2str(d(1)) f(i(3):i2(1)) 'mat'];
s=['RadAve_' site_name '_' file_date_str '.mat'];

eval(['save ' s ' Rclean;'])
%move to inital base directory
eval(['cd ' working_dir])  
      
%%% plot the result on a map of the array domain %%%
if CONST.goplot(1)==1
    figure(4); clf;
    HFR_DP_quickmap(ARRAY,working_dir);
    title(['Radave for Site:' RM.Site_name ', Time:' datestr(RM.time)])
    i=find(isnan(RADIAL.U(:,2))==0);
    m_quiver(RADIAL.LonLat(i,1),RADIAL.LonLat(i,2),RADIAL.U(i,2),RADIAL.V(i,2),10,'r')
%    m_quiver(RADIAL.LonLat(i,1),RADIAL.LonLat(i,2),RADIAL.U(i,1),RADIAL.V(i,1),10,'r')
    
    %%% if saving figures.
    if CONST.goprint(1)==1;
        eval(['cd ' outgoing_pics_file_dir])
        
        %save pictures of radials
        figure(4); hold on;
        t=title([s '      radave locations']); set(t,'interpreter','none');
        set(gcf,'paperposition',[.25 .25 5 5])
        i=find(s=='.');
        print('-djpeg','-r300',[s(1:i(1)-1) '_radave_' date '.jpg'])
        %move to inital base directory
        eval(['cd ' working_dir])
    end % if goprint
end% if goplot


%%%% convert existing, cleaned Rclean file to an ascii output in the
%%%% COS lluv format,

% %move to saving directory
% eval(['cd ' outgoing_radialavelluv_file_dir])
% 
[Rclean,HEAD]=HFR_DP_lera_RadAve2LLUV_v1(Rclean,CONST,patt,HEAD,SpecHead);
% 
% %move back to inital base directory
% eval(['cd ' working_dir])  


end %if RM is long enough to compute radial averages...
%%
end  % end jjj over fnames files


toc
diary off;
clear all; close all;

%%
return

