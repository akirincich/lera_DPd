%HFR_DP_master_LERA_XXXX.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HFR_DP_master_LERA_XXXX.m
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
%  See HFR_DP_LERA_SETUP_README.m for directions on how to set up a machine to
%  perform the processing
%
% v1        March 2017
% v2        December 2018
%
% Anthony Kirincich
% Woods Hole Oceanographic Institution
% akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%

clear
%%% load any paths for HFR_Progs, m_map, etc.  if needed...
%cd /home/hfrp/ ; HFR_DP_addpaths_hamlin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     begin user choices       %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set site/array information
CONST=[];

%%% the site name to process
CONST.site_name='XXXX';
CONST.Site_loc = [41+20.904/60 -(70+38.41/60)];
CONST.Site_bounds=[90 90+180];

%%% set the common constants for the HFR array this site is a part of:
ARRAY=[];
ARRAY.name='YYYY';
%%% array name, does not have to be four letters like radial sites, can be the
%%% same name for multiple radar sites that contribute to coverage within an area

%%% define the bounding box of the array, needed for any plotting operations.
ARRAY.min_lon=-(71+20/60);
ARRAY.max_lon=-(69+40/60);
ARRAY.min_lat=(40+40/60);
ARRAY.max_lat=(41+30/60);

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

%%% Where are the scripts, this is hard coded to be within the working directory
scripts_dir=[working_dir '/scripts/'];

%%% Where are the data  (see HFR_DP_SETUP_README.m for instructions
base_dir='/Users/anthony/Matlab/working/lera/lera_DP_testdata';
%base_dir='/Users/anthony/Data/LERA_process';
%base_dir='/Volumes/data/RadialSites';

%%% Where do the DP logs go, if they are enabled...
log_dir=['/Users/anthony/Matlab/working/LOGS'];
%log_dir=['/Volumes/data/LOGS'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% add needed paths
eval(['addpath ' scripts_dir]);
eval(['addpath ' working_dir]);
%%% move to inital base directory
eval(['cd ' working_dir])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set constants for the HFR
%%% how many solutions to process
CONST.Nmax=4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set parameters for spectral processing

%%% the final azimuthal resolution of the radial velocities, usually 5 deg less
CONST.bearing_width=1;  %in deg

%%% the minimum Signal to Noise ratio allowed (used in Radial averaging and ImageFOLs
CONST.snr_min=5;  %dB

%%%% the endpoints of the ave or RM pattern used for DF
CONST.pattern_ends=CONST.Site_bounds;  %These are equivalent and one should be removed

%%% for imageFOL, set user parameters
%%%  parameters are: [velD_change max_vel snr_min];  with vels defined in cm/s
%%%  see imageFOL for more details
%CONST.imageFOL_user_param=[20 100 CONST.snr_min];  %for a 25MHz site
%CONST.imageFOL_user_param=[40 300 CONST.snr_min];  %for a 5MHz site, viewing the Gulf Stream
CONST.imageFOL_user_param=[20 100 CONST.snr_min];  %for a 16MHz site
%CONST.imageFOL_user_param=[40 100 CONST.snr_min];  %or for a 16MHz site

%%% define radave RM thresholds... to be used to weed out bad radials before the radial average is made
%CONST.radave_thresholds=[snr_thresh angpeak_thresh angwidth_thresh(1) angwidth_thresh(2)];
%CONST.radave_thresholds=[CONST.snr_min 5 0 50];
CONST.radave_thresholds=[CONST.snr_min .001 0 200]; % for mle or wsf included

%%% define the peak threshold to use for the music doa peak finder,
%%% this is done in 10*log10(MS) space, making .05-.5 reasonable thresholds
%%% that scale with the values to allow small peaks to exist.
%CONST.doa_peak_thresh=.05;
%CONST.doa_peak_thresh=.25;
CONST.doa_peak_thresh=.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set processing choices

%%% Which first order line method to use
CONST.FOL_type='imageFOL';    %uses Image FOL methods
%CONST.FOL_type='codarFOL';    %imports the COS FOLs from the spectral file

%%% Which radial averaging method to output, both are estimated in the
%%% radave step, but only one can be sent out within the lluv file format.
CONST.radave_type='regular';    %follows arthmetic averaging of the radials

%  within the CONST.bearing_width area and error calcuation
%CONST.radave_type='metricQC';   %follows Kirincich et al 2012 to use thresholding
%  and power-weighted spatial means for rad ave

%%% which type of antenna manifold will be used
CONST.which_patt='Ideal';
%CONST.which_patt='Measured';

%%% which type of direction finding will be used
CONST.which_df='music';
%CONST.which_df='wsf';
%CONST.which_df='mle';

%%%% which type of method to determine the number of emitters...see Johnson and Dudgen, sec 7.3.5 for details
%CONST.which_Ns_meth='MDL';
%CONST.which_Ns_meth='AIC';
%CONST.which_Ns_meth='music_param';
CONST.which_Ns_meth='music_highest';

%%%% Sometimes the LERA Doppler spectra have lots of noise in them due to a
%%%%   nearby source of interference that maps into a small, and consist
%%%%   number of indices in each chirp.  We've found that masking and interpolating
%%%%   over small runs of the chirp returns gets rid of the bulk of the noise with
%%%%   little downside provide the run is 'small' (<10% of the chirp)
CONST.gohitch=1;    %find and remove a consistent 'hitch' in the chirps
%CONST.gohitch=0;    %process the spectra normally


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set processing constants, the only user selectable parts
CONST.PC=[];
CONST.PC.MaxRangeKM=100;        % the maximum ranges, in distance to be kept in spectra, from >0 to Max
%CONST.PC.SpecOverlapPct=50;     %for WOSA, style increased ensembles into CSE

CONST.PC.SpeclengthPnts=2048;   %defines the number of points in the spectra
CONST.PC.SpecOverlapPct=78;     %for WOSA, style increased ensembles into CSE
%CONST.PC.SpeclengthPnts=1024;   %defines the number of points in the spectra

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set plotting and printing switches
%plotting:
CONST.goplot=[1 0];   %where:
%  switch 1 controls the final plots of radial results, etc.  (1=plot, 0=do not plot)
%  switch 2 controls the intermediary plots within functions  (1=plot, 0=do not plot)
CONST.goprint=[0 0];
%  switch 1 prints the final plots of radial results, etc.  (1=print, 0=do not print)
%  switch 2 prints the intermediary plots within functions  (1=print, 0=do not print)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set how 'new' files to process are found

%%% Method 1: choose files based on the file time
CONST.files_to_process_method=1;
%%% if using this method, need to set the time bounds
%%% set the time frame to examine, only used by method 1
start_time=datenum(2018,9,1,0,0,0);
end_time=datenum(2018,11,1,0,0,0);
CONST.files_to_process_dates=[start_time end_time];

%%% Method 2: compare css files to RM files, and process new files
%CONST.files_to_process_method=2;


%%%%%%%%%%% end radial processing setup %%%%%%%%%%
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% running the processing scripts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% go and create those data folders that aren't already available
[site_name]=HFR_DP_LERA_createfilestructure_v1(CONST.site_name,base_dir);

% %%% setup logfile name for diary of the processing
% %%% this can be helpful if running in realtime
%  logfn=['LOG_v7_' mfilename '_'  datestr(now,30) '.log'];
%  diary([log_dir '/' logfn]);
%  tic

%%% for this site, prepare for processing
%%% sets up processing directories, loads site info and cal, and
%%% identifies files that need processing
HFR_DP_lera_ts2radial_prepwork_v4


%  %%% limit the number of files  on in any 1 pass of this processing.
%  %%% this can be helpful if running in realtime
%  if length(fnames)>5;
%      fnames=fnames(1:5);
%  end

%%
%%% now loop overeach file in fnames to process data
for jjj=1:1; %:length(fnames)
    
    %%% display the file to be processed
    filein=fnames{jjj}
    
    %%
    disp(filein)
    %%%% start processing steps list now that HEAD has been established.
    HEAD.ProcessingSteps={mfilename};
    %%% load the ts data for this file and make the ensemble-averaged spectra
    [SpecHead,data,mtime]=lera_time2spectra_v6(incoming_ts_file_dir,incoming_spectra_file_dir,filein,RC,HEAD,CONST);
    
    %%% save the spectra file
    %%%% get the output filename, based on the filein %%%%
    i=find(filein=='_');
    if strcmp(RC.radar_type,'MK2')==1;
        fileout=['CSE_' RC.SiteName '_' filein(i(2)+1:end)];
    elseif strcmp(RC.radar_type,'MK3')==1;
        fileout=['CSE_' RC.SiteName '_' filein(1:i(1)-1) '.mat'];
    end
    %%%%% save this data
    eval(['cd ' incoming_spectra_file_dir])
    save(fileout,'RC','SpecHead','data','mtime')
    eval(['cd ' working_dir])
    %%%%%%%%%%
 %%   
    if isempty(data)==1;  %something is wrong with this file...
        %save a placeholder RM file and move to the next one
        RM=[];
        %move to the RM directory
        eval(['cd ' outgoing_radialmat_file_dir])
        d=datevec(fdates(jjj));
        MM=num2str(d(2)); if d(2)<10; MM=['0' MM]; end;   dd=num2str(d(3)); if d(3)<10; dd=['0' dd]; end
        hh=num2str(d(4)); if d(4)<10; hh=['0' hh]; end;    mm=num2str(d(5)); if d(5)<10; mm=['0' mm]; end
        file_date_str=[num2str(d(1)) '_' MM  '_' dd '_' hh mm];
        
        s=['RM_' CONST.site_name '_' file_date_str '.mat'];
        note='This file had an issue and a blank RM file is being generated'
        eval(['save ' s ' RM RC SpecHead data mtime CONST note;'])
        %move to inital base directory
        eval(['cd ' base_dir])
        
        
    else  %%% continue as before ...
        
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
        %%%% run radial metric processing %%%
        HFR_DP_lera_spectra2radialmetric_process_v5    %update from 11/2018 for streamlined method selection
        
        %%%  again, output format of RM(jjj).data is similar to the radial metrics output
        %%%  with changes to the field contents to handle N=8
        %%%  mostly, only the results for one peak of one solution are given in
        %%%  each line,
        %%%
        %%%  cols   fields
        %%%   1-2    lon lat
        %%%   3-4    u v   (here nan as will not be used)
        %%%   5      flag    (here nan, used by COS but not here)
        %%%   6 range
        %%%   7 bearing
        %%%   8 vel
        %%%   9 direction
        %%%   10  rangecell
        %%%   11  dopcell
        %%%   12  angselect (which solution if isave out is given
        %%%   13  which peak of this solution is given
        %%%   14 musicpow(v)
        %%%   15 music_pkwidth
        %%%   16 musicDOASpeak
        %%%   17 spectraA3(snr)
        %%%   18-20 musicEigenRatio musicpowerRatio musicoffRatio
        %%%   21-27 which solutions were viable (1 yes 0 no, for 1-7 (need better description/name)
        %%%   28-35 Eigenval (1-8)
        %%%   36   DF_flag (see MUSIC algo above for explaination)
        
        
        %%% save this RM ouput...you don't have choice %%%
        %move to saving directory
        eval(['cd ' outgoing_radialmat_file_dir])
        d=datevec(fdates(jjj));
        MM=num2str(d(2)); if d(2)<10; MM=['0' MM]; end; dd=num2str(d(3)); if d(3)<10; dd=['0' dd]; end
        hh=num2str(d(4)); if d(4)<10; hh=['0' hh]; end; mm=num2str(d(5)); if d(5)<10; mm=['0' mm]; end
        file_date_str=[num2str(d(1)) '_' MM  '_' dd '_' hh mm];
        
        s=['RM_' CONST.site_name '_' file_date_str '.mat'];
        eval(['save ' s ' RM;'])
        %move to inital base directory
        eval(['cd ' base_dir])
        
        %%% plot the RM result on a map of the array domain %%%
        if CONST.goplot(1)==1
            figure(3); clf;
            HFR_DP_quickmap(ARRAY,working_dir);
            title(['RM for Site:' RM.Site_name ', Time:' datestr(RM.time)])
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
            s=['RadAve_' CONST.site_name '_' file_date_str '.mat'];
            
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
            %%%% IOOS required lluv format,
            
            % %move to saving directory
            eval(['cd ' outgoing_radialavelluv_file_dir])
            [Rclean,HEAD]=HFR_DP_lera_RadAve2LLUV_v1(Rclean,RC,CONST,patt,HEAD,SpecHead);
            % %move back to inital base directory
            eval(['cd ' working_dir])
            
            
        end %if RM is long enough to compute radial averages...
        %%
        
    end %if isempty(data)==1;  %if something is wrong with this file generate placeholder RM.
    
    
end  % end jjj over fnames files


toc
diary off;
clear all; close all;

%%
return

