function [RADIAL,HEAD]=HFR_DP_lera_RadialAverage_v2(RM,patt,HEAD,CONST);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [RADIAL,HEAD]=HFR_DP_RadialAverage_v5(RM,HEAD,CONST);
%
%  This script takes the Radial Metric output from HFR_DP and calculates
%   spatial averages over the user-defined azimuthal limits.  Note that
%   there is no temporal averaging available in the present package.
%
%  The averaged result is returned in a Matlab structure, named RADIAL,
%   that is a modified version of the HFR_Progs 'radial' structure to account
%   for the addition information returned with the radial metrics.
%
%  Two versions of the spatially averaged radials are reported:
%
%  (1) Arithmetic spatial mean
%
%  (2) Spatial means compute from QC'ed (threshold filtering) and MUSIC-derived
%      Power-weighted averaging following the work of Kirincich et al (2012)
%       and de Poalo et al (2015).   The threshold values are user-settable
%       and exist within the CONST.radave_thresholds set in the master
%       program.
%
%  NOTES:
%    (A) This script leans heavily on the HFR_Progs toolbox and recreates/modifies
%    some of the HFR_Progs scripts to handle the radial metrics output.  Thanks
%     Mike Cook and David Kaplan for their amazing initial efforts
%     developing the HFR_Progs toolbox and distributing it widely.
%
%    (B)  Background information on the file formats and the fields are
%    given in the script before the processing is done.
%
%  VERSIONS:
%     v1 create new for lera 8/2017
%     v2  7/2019 for gaussian averaged  and streamlined for readability
%
% Anthony Kirincich
% Woods Hole Oceanographic Institution
% akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%%%%%% background info on RM formats and what all the fields mean:  %%%%%%
%%% this is not critical for the script to run, but does give an interested
%%% user some background knowledge on what all the fields of the RM file
%%% are and how they might be able to be appled here.  Most of the info is
%%% found in the COS documention,
%
% The column output format of RM(jjj).data loosely follows COS RSv7 radial metric
%  file format with a few small changes to units and field content

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
% 17 spectrasnr(snr)
% 18-20 musicEigenRatio musicpowerRatio musicoffRatio
% 21-27 which solutions were viable (1 yes 0 no, for 1-7 (need better description/name)
% 28-35 Eigenval (1-8)
% 36   DF_flag (see above for explaination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Unpack processing thresholds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
snr_thresh=CONST.radave_thresholds(1)/5;  % 7/2019 change.
angpeak_thresh=CONST.radave_thresholds(2);
angwidth_thresh=CONST.radave_thresholds(3:4);

%%% define the list of ALL fields found within the data matrix coming from RM
vv={'LOND' 'LATD' 'VELU' 'VELV' 'VFLG' 'RNGE' 'BEAR' 'VELO' 'HEAD' 'SPRC' 'SPDC' ...
    'MSEL' 'MPKN' 'MPWR' 'DOAW' 'DOAP' 'MSNR' 'MEGR' 'MPKR' 'MDFR' ...
    'ISY1' 'ISY2' 'ISY3' 'ISY4' 'ISY5' 'ISY6' 'ISY7' ...
    'MEI1' 'MEI2' 'MEI3' 'MEI4' 'MEI5' 'MEI6' 'MEI7' 'MEI8' 'DFFG'};
%%% Basic List of 'radial' columns
cc = { 'LOND', 'LATD', 'RNGE', 'BEAR', 'HEAD', 'VELO' };
%%% full list of the non-radial fields that come out of RM
mm={'SPRC' 'SPDC' 'MSEL' 'MPKN' 'MPWR' 'DOAW' 'DOAP' 'MSNR' 'MEGR' 'MPKR' 'MDFR' 'ISY1' 'ISY2' 'ISY3' 'ISY4' 'ISY5' 'ISY6' 'ISY7' 'MEI1' 'MEI2' 'MEI3' 'MEI4' 'MEI5' 'MEI6' 'MEI7' 'MEI8' 'DFFG'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Other parameters to set, unlikely to change
%%% Set the Number of different estimates of the means to make
%nvel_est=4;
nvel_est=2;  %only return the power weighted average and straight mean

%%% set the file index to process
i5=1;  % I'm only going to do this for one file at a time, although you could
%  recast this function to do it for an array of RM structures with
%  many files. See the way Kaplan treats some of the radial
%  functions in HFR_Progs for an example.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  New in v2
% make the gaussian function for use below, makes a gaussian function with
% a 5 deg spread
ig=[-length(patt.Bear): median(patt.Bear-floor(patt.Bear)): length(patt.Bear) ] ;  %accounts for variance in the pattern decimal
Gaussian_alpha=160*1;
g=gausswin(length(ig),Gaussian_alpha)';
 i=find(g<.05); g(i)=nan;
if strcmp(CONST.radave_type,'gaus')==1;
    disp('doing gaussian weighted averaging')
if CONST.goplot(2)==1;
    figure(10); clf; plot(ig,g); hg; length(find(g>.5))/2
    title(['using Gaussian with size ' num2str(length(find(g>.5))/2)])
end
else
    disp('doing power weighted averaging')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% place the RM data into a modified RADIAL structure
R = RADIALstructv7;   %get a blank structure.
R.FileName = RM(i5).fname;
%print to screen so we know what is happening
fprintf('Processing radial averages for RM of: %s\n',R.FileName)
% Give TimeStamp and matvars same size as FileName for consistency just
% in case we return early.
R.TimeStamp = RM(i5).time;   %time
R.TimeZone = 'GMT';          %timezone

% Add processing step to the structure
R.ProcessingSteps{end+1} = [RM(i5).fname ' to RM with HFR_spectra2radialmetrics_process_v?'];
R.SiteOrigin = RM(i5).Site_loc([2 1]);   % Origin Lat,Lon,
R.SiteName = RM(i5).Site_name;           % Site name
R.Type = patt.Pattern_type;                         % PatternType

%%% find where the data of interest is:
%%%   This might seem overly complex now, but it will allow for the RM file
%%%   to be flexible in format without breaking the script.
%%% Get the indices of critical fields for radial info only
II = [];
for k = 1:length(cc);
    ii = strmatch( cc{k}, vv, 'exact' );
    if isempty(ii)
        feval( warnfunc, [ mfilename ':MISSING_ESSENTIAL_DATA_COLUMN' ], ...
            [ 'One or more of the required data columns cannot be found.' ...
            ' Returning empty structure.' ] );
        %  return
    else;  II(k) = ii;
    end
end

% Get pieces of data we want if they are not empty.
R.LonLat = RM(i5).data( :, II(1:2) );
R.RangeBearHead = RM(i5).data( :, II(3:5) );
R.RadComp = RM(i5).data(:,II(6))*100;  %convert to cm/s for rest of processing stream;

% !!!!!!!!!!!!!!!!!!!!!!!!!!
% Change direction to cartesian convention.  Note, however, that bearing
% will still point away from radar and heading will still point towards radar.
R.RangeBearHead(:,2:3) = true2math( R.RangeBearHead(:,2:3) );
% !!!!!!!!!!!!!!!!!!!!!!!!!!

%%% Get the indices of critical fields for radial metrics
JJ = [];
for k = 1:length(mm);
    ii = strmatch( mm{k}, vv, 'exact' );
    if isempty(ii)
        feval( warnfunc, [ mfilename ':MISSING_ESSENTIAL_DATA_COLUMN' ], ...
            [ 'One or more of the required data columns cannot be found.' ...
            ' Returning empty structure.' ] );
        %  return
    else; JJ(k) = ii;
    end
end

%output the whole thing as a large matrix
R.Metric=RM(i5).data(:,JJ);
% Add currently unused Flag (1 is for original radials).
R.Flag = ones( size(R.RadComp) );

%%

%%%%%%% new parts that convert the radial metric data to regular radials %%%%% %%%%%
if isempty(R.Metric)==0 & isempty(find(~isnan(R.Metric(:))==1))==0  %%% skip files that are empty
    %%
    %%%% doing this range bin separately for each file will allow variable
    %%%% ranges to be combined into the same timeseries...although later
    %%%% averaging will have to account for this as well.
    
    %%% make aver bin range points
    [range]=unique(R.RangeBearHead(:,1));
    %for the Bear and range...make to latlon using m_fdist
    FLonLat=[]; FRangeBearHead=[];
    a=1;
    for j=1:length(patt.BearT)
        for i=1:length(range)
            [FLonLat(a,1),FLonLat(a,2),FLonLat(a,3)]=m_fdist(R.SiteOrigin(1),R.SiteOrigin(2),patt.BearT(j),range(i).*1000);
            %also put out the current range and Bear
            FRangeBearHead(a,:)=[range(i) patt.BearT(j) FLonLat(a,3)];
            a=a+1;
        end
    end
    %%% sometimes the longitude comes out greater than +/180deg, im not
    %%% sure why.  and its not wrong mbut interfers with the plotting, fix
    i=find( FLonLat(:,1) >180);    FLonLat(i,1)=FLonLat(i,1)-360;
    FRangeBearHead(:,2)=true2math(FRangeBearHead(:,2));
    FRangeBearHead(:,3)=true2math(FRangeBearHead(:,3));
    
    %  %Test look at this file layout
    if CONST.goplot(2)==1
        figure(10); clf;
        plot(FLonLat(:,1),FLonLat(:,2),'r.');
        title([R.FileName ' Radave locations ';],'interpreter','none')
    end
    
    %%% do same for the bearing midpoints at the largest range.
    i=find(range==max(range));
    FmLonLat=[];
    a=1;
    for j=1:length(patt.BearT_ends)
%        [FmLonLat(a,1),FmLonLat(a,2),FmLonLat(a,3)]=m_fdist(R.SiteOrigin(1),R.SiteOrigin(2),patt.BearT_ends(j),range(i).*1000);
        [FmLonLat(a,1),FmLonLat(a,2),FmLonLat(a,3)]=m_fdist(R.SiteOrigin(1),R.SiteOrigin(2),patt.BearT(j),range(i).*1000);
        a=a+1;
    end
    %%% sometimes the longitude comes out greater than +/180deg, im not
    %%% sure why.  and its not wrong mbut interfers with the plotting, fix
    i=find( FmLonLat(:,1) >180);    FmLonLat(i,1)=FmLonLat(i,1)-360;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Start outgoing radial structure
    RADIAL(i5)=R;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% new section to process this data for the given patt file
    % DOA Peak Value                               
    % DOA Function Width (at half power)            
    % Signal Power (a measure of the eigenvalue)    
    % Antenna 1,2,3 SNR for this doppler cell       
    
    %Build output structure with gridded products....
    RADIAL(i5).LonLat=FLonLat(:,1:2);
    RADIAL(i5).RangeBearHead=FRangeBearHead;
    RADIAL(i5).RadComp=nan*ones(length(FLonLat),nvel_est); RADIAL(i5).Error=RADIAL(i5).RadComp;
    RADIAL(i5).U=RADIAL(i5).RadComp; RADIAL(i5).V=RADIAL(i5).RadComp;
    RADIAL(i5).Flag=nan*ones(length(FLonLat),2);
    RADIAL(i5).Metric=nan*ones(length(RADIAL(i5).LonLat),6);   %for carrying scripps-type average metrics forwards
    RADIAL(i5).ProcessingSteps(length(RADIAL(i5).ProcessingSteps)+1)={'weighted average in HFRRMetricLoad_process'};
    
    %%
    for jj=1:length(range);     %loop over each range cell,
    %%
    ii=find(R.Metric(:,1)==jj);   %find all data that falls in this range 
         % in lera files, the bearing selected is given outright in
         % each line
         mbearing=R.RangeBearHead(ii,2);
        
   %     RAM=[]; % clear a structure that would hold all the individal groups  %This is not needed in v2
        % of radial for this range circle
        %%% for each ave bearing window, in this range cell, place all the found velocities within a new matrix
        for kk=1:length(patt.Bear);
              
            %%% all the points that match this bearing exactly
           i3=find(mbearing<max(patt.Bear_ends(:,kk)) & mbearing>=min(patt.Bear_ends(:,kk))) ; 
            %find the relbear of all points on this range circle
               %%% for the gaussian smoother
               relbear=mbearing-patt.Bear(kk);
               %get rid of bearings far out of domain
               i=find(abs(relbear)<length(patt.Bear));  i1=ii(i); relbear=relbear(i);
               %get the indices of the ig where data is found
               irelbear=[];
               for j=1:length(relbear)
                   irelbear(j)=find(ig==relbear(j));
               end      

                %%%% find the azimuthal deg grid point that matches this bin
                i2=find(RADIAL(i5).RangeBearHead(:,1)==range(jj));
                i22=find(RADIAL(i5).RangeBearHead(i2,2)==patt.Bear(kk));
                %so use i2(i22) to place the data in the outgoing structure

                %setup with blank returns
                a=[nan nan]; b=[nan nan]; met=nan.*ones(1,6); ave_solution=nan;
     if ~isempty(i3)==1;  %if data exists at this bearing
     %%%%% for the data at this bearing
                %get power, volt, ang_width, and ang_peak for each velocity in list
                sp=R.RadComp(ii(i3));
                %%% in the RM files 10:12 are the signal voltages, not powers
                volt=R.Metric(ii(i3),5); %RAM(kk).Metric(:,5);
                power=20*log10(abs(volt)) + 200;   % The db has to be pos, for power weighting     
                ang_width=R.Metric(ii(i3),6); %RAM(kk).Metric(:,6);
                ang_peak=R.Metric(ii(i3),7); %RAM(kk).Metric(:,7); % what is the real unit of this number?
                ant_snr=R.Metric(ii(i3),8); %RAM(kk).Metric(:,8);
                
    %%%%%%% time to cut data with bad... for regular data
                isnr=find(ant_snr < snr_thresh);   
                ant_snr(isnr,:)=nan;
                isnr_a=find(isnan(mean(ant_snr,2))==0);
                isnr_b=find(isnan(mean(ant_snr,2))==1);  %for flag               
                %%%ang_width
                iaw_a=find(ang_width > angwidth_thresh(1) & ang_width < angwidth_thresh(2));
                iaw_b=find(ang_width < angwidth_thresh(1) | ang_width > angwidth_thresh(2)); %for flag
                %%% ang_peak
                iap_a=find(ang_peak > angpeak_thresh);
                iap_b=find(ang_peak < angpeak_thresh);  %for flags                
                %%%%%% find the radials that pass all the tests
                igood=intersect(isnr_a,intersect(iaw_a,iap_a));
               
                %%% make flags for this
                %RADIAL(i4).Flag(i2(i22),2)=[length(igood)*1e6 + length(isnr_b)*1e4 + length(iaw_b)*1e2+length(iap_b)];
                %%%  [ #of_ang_peak_cut (2digits)  #of_ang_width_cut (2digits)    #of_snr3_cut  # of igood (2digits)  ];
                RADIAL(i5).Flag(i2(i22),2)=[ length(iap_b)*1e6 + length(iaw_b)*1e4 + length(isnr_b)*1e2 + length(igood) ];
                
                   if length(igood)>=1   %continue to make estimates of the average radial velocity                       
                        %%establish data and weights arrays
                        d=sp(igood);
                       % w1=volt(igood)./sum(volt(igood));   %voltage weighting
                        w1=power(igood)./sum(power(igood));   %power weighting
                        %average
                        a=[nanmean(d) sum(d.*w1)];                       
                        if length(igood)>=3
                            %do std devation   %%%% see http://en.wikipedia.org/wiki/Weighted_mean
                            %%% w1 is normalized to sum(w1)=1, thus an 'unbiased' weighted variance
                            %%%      is not possible and sig^2 = sum(w1.*(d-a(2).^2)./V1  where V1 =sum(w1)=1
                            b= sqrt( [sum( (d-a(1)).^2)./(length(igood)-1) ...
                                sum( w1.*((d-a(2)).^2))./1 ] );
                        end
                   %return other average metrics for later QA/QC
                    % (1) DOA Peak Value  (2) DOA Function Width (at half power) (3) Signal Power (4-6) ant1-3 snr
                    met=[mean(ang_peak(igood)) mean(ang_width(igood)) mean(10*log10(nanmean(volt(igood)))+ 60) nan nan mean(ant_snr(igood),1) ];
                    ave_solution=nanmean(R.Metric(ii(i3(igood)),3));
                    end
                 
     end % if ~isempty(i3)==1;  %if data exists at this bearing
                      
                    
 if strcmp(CONST.radave_type,'gaus')==1 & ~isempty(irelbear)==1;  %if some data exists along this range cell proceed
%%%%% for the data within this relative gaussian window range
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %%% for the data along the range circle, find the properties to ave
                    %%% with a gaussian
                    dall=R.RadComp(i1);
                    powerall=20*log10(abs(R.Metric(i1,5))) + 200 ;   % The db has to be pos, this is done in the next step for sure,
                    %powerall=powerall-(min(powerall)*[1-sign(min(powerall))*.1]); 
                    powerall=powerall./max(powerall);  %normalize weight
                    powerall(powerall<0)=nan;                    %remove low values 
                    dpeak=10*log10(R.Metric(i1,7));
                    ant_snrall=R.Metric(i1,8);
                    %gaussian weight
                    gv=powerall.*g(irelbear)';
         %           igoodall=find(~isnan(dall)==1 & dpeak > CONST.doa_peak_thresh & isnan(gv)==0 & ant_snrall > snr_thresh);
                    igoodall=find(~isnan(dall)==1 & isnan(gv)==0 & ant_snrall > snr_thresh);
                    %normalize by data that is available. 
                    gv=gv(igoodall)./nansum(gv(igoodall));                     

                    if isempty(find(gv>1 | gv<0))==0;
                        disp('there is something wrong with the weight function')
                    end
              
                 if length(igoodall)>=1   %continue to make estimates of the average radial velocity
                       %average
                        a(2)=[sum(dall(igoodall).*gv)];                        
                        if length(igoodall)>=3  %how the std dev is estimated
                            %do std devation   %%%% see http://en.wikipedia.org/wiki/Weighted_mean
                            %%% w1 is normalized to sum(w1)=1, thus an 'unbiased' weighted variance
                            %%%      is not possible and sig^2 = sum(w1.*(d-a(2).^2)./V1  where V1 =sum(w1)=1
                            b(2)= sqrt( sum( gv.*((dall(igoodall)-a(2)).^2))./1  );
                        end  %igood length >3 or not  
                       %return other average metrics for later QA/QC
                    % (1) DOA Peak Value  (2) DOA Function Width (at half power) (3) Signal Power (4-6) ant1-3 snr
                    met=[mean(dpeak(igoodall)) nan mean(powerall(igoodall)) nan nan mean(ant_snrall(igoodall),1) ];
                 ave_solution=nanmean(R.Metric(i1(igoodall),3));
                 end
             % pause   
        end  %if go irelbear
        
     %   if abs(diff(a))>50; andfa; end
        
                    RADIAL(i5).RadComp(i2(i22),:)=a;
                    RADIAL(i5).Error(i2(i22),:)=b;
                    %whats the average solution used here 1,2, or 3?
                    RADIAL(i5).Flag(i2(i22),1)=ave_solution; %nanmean(R.Metric(igood,3)); %nanmean(RAM(kk).Metric(igood,3));
                    %return other average metrics for later QA/QC
                    % (1) DOA Peak Value  (2) DOA Function Width (at half power) (3) Signal Power (4-6) ant1-3 snr
                    RADIAL(i5).Metric(i2(i22),:)=met; %[mean(ang_peak(igood)) mean(ang_width(igood)) mean(10*log10(nanmean(volt(igood)))+ (-40 - 5.8)) nan nan mean(ant_snr(igood),1) ];
                    
                    
                    %%%% plot results, if we are in the weeds with looking
                    %%%% at all steps of the processing
                    if length(i3)>5 & CONST.goplot(2)==1;
                        
                        angle=mbearing(i3);
                        ant3snr=ant_snr;
                        %regroup by angle
                        [s si]=sort(angle);
                         
                        figure(1); clf
                        subplot(311);
                        plot(1:length(angle),angle(si),'ko','markerfacecolor','k'); hold on;
                        title('angle'); ylabel('deg'); %xlabel('# radials')
                        axis([0 length(angle)+1 min(angle)-1 max(angle)+1])
                        
                        
                        s1=subplot(312);
                        set(s1,'box','off','xlim',[0 length(angle)+1],'xtick',[0:length(angle)],...
                            'ylim',[0 max(ant3snr(si))])
                        l1=line(1:length(angle),ant3snr(si),'color','k','parent',s1); hold on
                        l2=line(1:length(angle),ant3snr(si),'color','k','marker','.','parent',s1);
                        title('ant3snr and speed')
                        ylabel('ant3snr');
                                                
                        a1=axes('position',get(s1,'position')); %speed
                        set(a1,'yaxislocation','right','xaxislocation','top','box','off',...
                            'color','none','xlim',[0 length(angle)+1],'xtick',[0:length(angle)],...
                            'ylim',[-35 35],'xticklabel',[])
                        l2=line(1:length(angle),sp(si),'color',[.5 .5 .5],'parent',a1); hold on;
                        l2=line(1:length(angle),sp(si),'color',[.5 .5 .5],'marker','.','parent',a1); hold on;
                        pp=line([0 length(angle)+1],[a(1) a(1)],'color','k','parent',a1(1));
                        
                        pp=line([0 length(angle)+1],[a(2) a(2)],'color','k','linestyle','--','parent',a1);
                        ylabel('speed cm/s');
                                                
                        s1=subplot(313);
                        set(s1,'box','off','xlim',[0 length(angle)+1],'xtick',[0:length(angle)],...
                            'ylim',[-50 0])
                        l1=line(1:length(angle),power(si),'color','k','parent',s1); hold on
                        l2=line(1:length(angle),power(si),'color','k','marker','.','parent',s1);
                        title('power and speed')
                        ylabel('power dB');
                        xlabel('# radials');
                        
                        a1=axes('position',get(s1,'position')); %speed
                        set(a1,'yaxislocation','right','xaxislocation','top','box','off',...
                            'color','none','xlim',[0 length(angle)+1],'xtick',[0:length(angle)],...
                            'ylim',[-35 35],'xticklabel',[])
                        l2=line(1:length(angle),sp(si),'color',[.5 .5 .5],'parent',a1); hold on;
                        l2=line(1:length(angle),sp(si),'color',[.5 .5 .5],'marker','.','parent',a1); hold on;
                        pp=line([0 length(angle)+1],[a(1) a(1)],'color','k','parent',a1(1));
                        
                        pp=line([0 length(angle)+1],[a(2) a(2)],'color','k','linestyle','--','parent',a1);
                        ylabel('speed cm/s');
                        
                        pause  %if your plotting this, you wish to look at every single result
                    end  % if length(i3)  plot
               % end  %if length(igood)==                
            %end %isempty(i)
        end %for kk
        
    end %jj length range
    
    %%
    
    %%%% make U,V velocities for plotting %%%%
    H=math2true(RADIAL(i5).RangeBearHead(:,3));
    RADIAL(i5).U(:,:)=RADIAL(i5).RadComp(:,:).*sin(H*ones(1,nvel_est).*pi/180);
    RADIAL(i5).V(:,:)=RADIAL(i5).RadComp(:,:).*cos(H*ones(1,nvel_est).*pi/180);
    
    HEAD.ProcessingSteps{end+1}=mfilename;
    RADIAL.ProcessingSteps=HEAD.ProcessingSteps;
    
end  %isempty(R.Metric)==0

