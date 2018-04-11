%HFR_DP_lera_spectra2radialmetrics_process_v1.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HFR_DP_lera_spectra2radialmetrics_process_v1.m
%
%  This script runs the bulk of the processing steps used to derive radial 
%  velocity estimates from the spectra observed from lera-type HF radars
%  
%  In short, this script performs the following operations using the
%  user-defined choices and constants described by HFR_DP_master_XXXX.m,
%  including:
%
%  -Set up the outgoing radial metrics file for this spectral estimate
%  -Form the file-specific meta data critical
%      to estimating radial velocities.
%  -Establish the first order limits (FOLs)
%  -Use MUSIC to perform direction finding on the FO data 
%  -Finallize the output MATLAB Structure, that contains the 'radial
%      metric' results 
%
% INPUT:  As this is a script, not a function, everything in the workspace
%           is passed to the script for its use
%
% OUTPUT: As this is a script, not a function, everything is passed back,
%          but the key addition is:
%
%   RM  --  A matlab structure that matrix of radial metric output,
%           including the following fields:
% 
% RM.Site_name            The 4 character site name
% RM.Site_loc             [Lon Lat] of the radar
% RM.fname                The file name processed
% RM.time                 Matlab time of the observed spectra/radials
% RM.rang                 Ranges (km)
% RM.c_vel                Doppler velocity of the spectral estimate (m/s)           
% RM.Stats                the fraction of single and duel angle solutions reported
% RM.data                 *** See Below ***
% RM.FOregi               the spectral indices that are processed
% RM.ProcessingSteps      A cell array list of all the scripts/functions
%                           used to get to this point in the processing
%
% Output format of RM.data is significantly different than 
% codar-type output radial metric files
% cols   fields
% 1-2    lat lon
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
% 36   DF_flag (see above for explaination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Version: 
% v1       August 2017
%
%
% Anthony Kirincich
%   WHOI-PO
%   akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%% mark that we are using this file to process...
HEAD.ProcessingSteps{end+1}=mfilename;

    %create output structure
    RM=[];
    RM.Site_name=HEAD.SiteName;
    RM.Site_loc=HEAD.Site_loc;
    RM.fname=[];
    RM.time=[];
    RM.rang=[];
    RM.c_vel=[];
    RM.Stats=[];
    RM.data=nan.*ones(1,36);
    RM.FOregi=[];

%%%%%%%%%%%%%%%%%%%%%%
%do FOLs
%%%%%%%%%%%%%%%%%%%%%%

%%% set last constant
CONST.v_incr=median(diff(SpecHead.c_velc));
CONST.whichant_FOL=8; % can be a single antenna or a group of antennas
%%% run imageFOLS
%[FOreg,FOregi,Alims,HEAD,DN_out]=lera_imageFOLs_v1(data.a33,SpecHead.iFBraggc,CONST.v_incr,CONST.imageFOL_user_param,CONST.goplot,HEAD);
[FOreg,FOregi,Alims,HEAD,DN_out]=lera_imageFOLs_v2(data,CONST.whichant_FOL,SpecHead.iFBraggc,CONST.v_incr,CONST.imageFOL_user_param,CONST.goplot,HEAD);

%%% test this function
%%% ant3_volt=data.a33;
% iant=CONST.whichant_FOL;
% iFBragg=SpecHead.iFBraggc;
% v_incr=CONST.v_incr;
% user_param=CONST.imageFOL_user_param;
% goplot=[1 1]; %CONST.goplot;

if isempty(FOreg)==1;
    disp(['error...incorrect FOL_type: ' char(CONST.FOL_type) ' not identified'])
    asdfasdfsdf
end

%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% process FO data with music if the file has enough FO data to be
%%% interesting
if length(FOregi)>100
    tic ; tag=0;
%    if CONST.goplot(2)==1; CONST.goplot(2)=0; tag=1; end %turn off extreme plotting for this fn
    [R, HEAD]=HFR_DP_lera_music_v1(data,FOregi,HEAD,SpecHead,patt,CONST);
%    if tag==1; CONST.goplot(2)=1; end  %turn it back on
    toc
else
    R=nan.*ones(1,31);
end
        
%%
%only move forward if have a significant number of returns...i.e.
%length(R) > 500
if length(R)>500;
    
    %%% Finish radial metric file creation,
    %%%  calc lon lat u v flag for each line
    r=nan*ones(length(R),5);
    
    x=R(:,1)*[1 1].*[sind(R(:,2)) cosd(R(:,2))];   %reverse sin cos to account for true vs math coord
    [lon,lat]=km2lonlat([HEAD.Site_loc(2).*ones(length(R),1)],[HEAD.Site_loc(1).*ones(length(R),1)],x(:,1),x(:,2));
    r(:,1:2)=[lon lat];
    %find U,V
    r(:,3:4) = [R(:,3).*cosd(true2math(R(:,4)))    R(:,3).*sind(true2math(R(:,4)))];
    
    %         %%%% test range, bearing estimates of x,y in km
    %         d=[45 135 225 315];
    %         dm=true2math(d)
    %         x=R(1,1).*[sind(d)' cosd(d)']  %reverse sin cos to account for true vs math coord
    %         xm=R(1,1).*[cosd(dm)' sind(dm)']
    
    %%
    
    %place results into metrics file
    RM.fname=char(fnames(jjj));
    RM.time=mtime+15/60/24;  %move time to the midpoint of the file
    RM.rang=SpecHead.RangeKm;
    RM.c_vel=SpecHead.c_velc;
    RM.FOregi=FOregi;
    i=find(R(:,7)==1);
    
    % stats on the fraction of single and duel angle solutions being
    % reported
    for i=1:CONST.Nmax
        RM.Stats(i)=[length(find(R(:,7)==i))./length(R)];
    end
    
    RM.data=[r R];
    
    RM(1).ProcessingSteps=HEAD.ProcessingSteps;
%    RM.patt_UUID=patt.UUID;
    
    %again, output format of RM(jjj).data is similar to COS RSv7 radial metric files
    %with changes to the field contents to handle N=8
    %mostly, only the results for one peak of one solution are given in
    %each line,
% cols   fields
% 1-2    lat lon
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

    
else  %if there wasn't enough data to make a real file
    RM(1).fname=char(fnames(jjj));
    RM(1).time=mtime+15/24/60;
    RM(1).rang=SpecHead.RangeKm;
    RM(1).c_vel=SpecHead.c_velc;
    %   i=find(R(:,7)==1);
    %   RM(jjj).Stats=[length(i)./length(R) (length(R)-length(i))./length(R)];
    %   RM(jjj).data=[r R];
    RM(1).Stats=nan;
    RM(1).data=nan.*ones(1,36);
    RM(1).ProcessingSteps=HEAD.ProcessingSteps;
    RM(1).patt_UUID=patt.UUID;
    
end  %length(R) >500

%%% go back main script to save RM and possibly the figures




