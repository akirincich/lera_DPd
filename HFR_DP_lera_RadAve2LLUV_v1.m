function  [Rclean,HEAD]=HFR_DP_lera_RadAve2LLUV_v1(Rclean,RC,CONST,patt,HEAD,SpecHead);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HFR_DP_radave2lluv.m
%
% Script uses PATT, SpecHead, HEAD, and Rclean to fill out the required 
%  elements of the COS lluv file format, writes a new file to the current 
%  directory with the outfilename:
%
%   ['RDLm' '_' char(Rclean.SiteName) '_' tt '.ruv'];
%
%  where tt is the time of the data.  The file is written and saved inside
%  of the function into the current directory.  This script doesn't alter 
%  Rclean or HEAD in any way, just returns them for something to return.
%
%  Version:
%  -v1  adapted from cos V5 files but made sufficent with lera DF processing  
%
%  Anthony Kirincich
%  WHOI-PO
%  akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
Rclean.ProcessingSteps{end+1}=mfilename;

%new constants for output to lluv
PROC_CONST=[];
PROC_CONST.MergedCount=1;
PROC_CONST.RadialMinimumMergePoints=1;
PROC_CONST.Temp_qual=999;
PROC_CONST.VelocityMaximum=999;
PROC_CONST.VelocityMimimum=999 ;
PROC_CONST.TemporalCount=1;
PROC_CONST.COS_Flag=999;

if strcmp(CONST.radave_type,'regular')==1;
    PROC_CONST.iwhich_ave=1;
elseif strcmp(CONST.radave_type,'metricQC')==1;
    PROC_CONST.iwhich_ave=2;
end

%%% use Rclean.U to figure out which rows are good.
ikeep=find(isnan(Rclean.RadComp(:,PROC_CONST.iwhich_ave))==0);


% %%%%%%%% setup proper outputfile name  and directory %%%%%%%%%%%
%make proper file name 
if strcmp(CONST.which_patt,'meas')==1
    f_type='RDLm';
elseif strcmp(CONST.which_patt,'ideal')==1
    f_type='RDLi';
end
%set proper time stamp    
tt=datestr(Rclean.TimeStamp,31); tt=tt(1:16);
i=find(tt=='-' | tt==' '); tt(i)='_';
i=find(tt~=':'); tt=tt(i);
outfilename=[f_type '_' char(Rclean.SiteName) '_' tt '.ruv'];

%%% lets go!
fid=fopen(outfilename,'w');
fprintf(fid,'%%CTF: 1.00\n%%FileType: LLUV rdls "RadialMap"\n%%LLUVSpec: 1.0 2018 11 28\n%%Manufacturer: UH LERA\n');   %new for lera

fprintf(fid,'%%Site: %s ""\n',char(Rclean.SiteName));
tt=datestr(Rclean.TimeStamp,31); i=find(tt=='-' | tt==':'); tt(i)=' ';
fprintf(fid,'%%TimeStamp: %s\n',tt);

%i=find(Rclean.TimeZone==' ');
%tzone=Rclean.TimeZone(1:i(1)-1);
%toff=Rclean.TimeZone(i(1)+1:i(2)-1);
%there is never daylight savings time
%fprintf(fid,'%%TimeZone: "%s" +%s 0\n',tzone,toff);  %it should always be GMT.
fprintf(fid,'%%TimeZone: "UTC" +0 0\n');  %it should always be UTC.  (GMT has DST)
fprintf(fid,'%%TimeCoverage: %2.4f Minutes\n',SpecHead.nchirps*SpecHead.Tr/60);
fprintf(fid,'%%Origin: %3.7f %3.7f\n', CONST.Site_loc);
fprintf(fid,'%%GreatCircle: "WGS84" 6378137.000  298.257223562997''\n');
fprintf(fid,'%%GeodVersion: "CGEO" 1.57  2009 03 10\n');
fprintf(fid,'%%LLUVTrustData: all %%%% all lluv xyuv rbvd\n');

fprintf(fid,'%%RangeStart: 1\n');
fprintf(fid,'%%RangeEnd: %2.0f\n',length(SpecHead.RangeKm));
%%% COS format for the LLUV file is to report the results in meters if less
%%% then 500 m, but KMeters if over
if SpecHead.RangeResKm < 0.5
    fprintf(fid,'%%RangeResolutionMeters: %1.7f\n',SpecHead.RangeResKm.*1000);
elseif SpecHead.RangeResKm > 0.5
    fprintf(fid,'%%RangeResolutionKMeters: %1.7f\n',SpecHead.RangeResKm);
end


fprintf(fid,'%%RxAntConfig: %s\n',RC.RxAntConfig);  %new to lera
fprintf(fid,'%%TxAntConfig: %s\n',RC.TxAntConfig);  %new to lera
fprintf(fid,'%%RxAntennaBearing: %3.1f True\n',patt.Array_bearing);   %new to lera
fprintf(fid,'%%TxAntennaBearing: %3.1f True\n',RC.Tx_bearing);  %new to lera

fprintf(fid,'%%ReferenceBearing: 0 True\n');
fprintf(fid,'%%AngularResolution: %1.0f Deg\n',median(diff(patt.Bear)));
fprintf(fid,'%%SpatialResolution: %1.0f Deg\n',median(diff(patt.Bear)));
fprintf(fid,'%%patternType: %s\n',patt.Pattern_type);

tt=datestr(patt.date,31); i=find(tt=='-' | tt==':'); tt(i)=' ';
fprintf(fid,'%%PatternDate: %s\n',tt);

fprintf(fid,'%%PatternResolution: %1.1f deg\n',median(diff(patt.Bear)));
fprintf(fid,'%%PatternSmoothing: %2.1f deg\n',median(diff(patt.Bear)));
%fprintf(fid,'%%PatternUUID: %s\n',PATT.UUID);
%fprintf(fid,'%%PatternUUID: n/a\n');
    

%if SpecHead.bSweepUp==1
fprintf(fid,'%%TransmitCenterFreqMHz: %2.6f\n',SpecHead.Fc/1e6);  %only Fc is needed for lera
%else 
%    fprintf(fid,'%%TransmitCenterFreqMHz: %2.6f\n',SpecHead.fStartFreqMHz + -1.*SpecHead.fBandwidthKHz/1000/2);
%end
fprintf(fid,'%%DopplerResolutionHzPerBin: %1.9f\n',SpecHead.delta_f);
%fprintf(fid,'%%BraggSmoothingPoints: %2.0f\n',PROC_CONST.BraggSmoothingPoints);
fprintf(fid,'%%CurrentVelocityLimit: %3.1f\n',CONST.imageFOL_user_param(2));

%fprintf(fid,'%%BraggHasSecondOrder: %2.2f\n',PROC_CONST.BraggHasSecondOrder);
%fprintf(fid,'%%RadialBraggPeakDropOff: %2.2f\n',PROC_CONST.RadialBraggPeakDropOff);
%fprintf(fid,'%%RadialBraggPeakNull: %2.2f\n',PROC_CONST.RadialBraggPeakNull);
%fprintf(fid,'%%RadialBraggNoiseThreshold: %2.2f\n',PROC_CONST.RadialBraggNoiseThreshold);

%%% use the max real part of each antenna as the ampfactors, this would
%%% change from the ideal pattern of 1 to something else if the pattern was
%%% non-ideal 
ampfactors=max(abs(real(patt.A))');
fprintf(fid,'%%PatternAmplitudeCorrections: %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n',ampfactors);
%%% use the imag part to get the phase offsets
%%% as the reference is always 0 find the minimum phase value
phasefactors=min(abs(imag(patt.A))');
fprintf(fid,'%%PatternPhaseCorrections: %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n',phasefactors);

%fprintf(fid,'%%PatternAmplitudeCalculations: %1.4f %1.4f\n',HEAD.Ampfacs_ideal1_2_meas1_2(3:4));
%fprintf(fid,'%%PatternPhaseCalculations: 999 999\n');
%fprintf(fid,'%%RadialMusicParameters: %2.3f %2.3f %2.3f\n',HEAD.Musicparams123_gsmoothwidth_deg_smearwidth_velthresh(1:3));

fprintf(fid,'%%MergedCount: 1\n');
fprintf(fid,'%%RadialMinimumMergePoints: %2.0f\n',PROC_CONST.RadialMinimumMergePoints);
%fprintf(fid,'%%MergeMethod: Arithmetic Mean or Weighted Mean\n');

%specific to the HFR_DP methods.
fprintf(fid,'%%RadialSpatialAveragingMethod_HFR_DP: %s\n',CONST.radave_type);
fprintf(fid,'%%FirstOrderCalcMethod: %s\n',CONST.FOL_type);
fprintf(fid,'%%FirstOrderCalcParameters: %2.3f %2.3f %2.3f\n',CONST.imageFOL_user_param);

fprintf(fid,'%%RadialVelocityEstimateMethod: DirectionFinding\n');
fprintf(fid,'%%DirectionFindingMethod_HFR_DP: %s\n',CONST.which_df);
fprintf(fid,'%%DirectionFindingEmitterDetection_HFR_DP: %s\n',CONST.which_Ns_meth);
  fprintf(fid,'%%RadialMusicDOAPeakThreshold: %2.3f\n',CONST.doa_peak_thresh);
if strcmp(CONST.which_Ns_meth,'music_param')==1
  fprintf(fid,'%%RadialMusicParameters: %2.3f %2.3f %2.3f\n',HEAD.Musicparams123);
end

fprintf(fid,'%%PatternMethod: 1 PatternVectors\n');
fprintf(fid,'%%TransmitSweepRateHz: %2.4f\n',1/SpecHead.Tr);
fprintf(fid,'%%TransmitBandwidthKHz: %3.6f\n',RC.BW/1000);
fprintf(fid,'%%SpectraRangeCells: %4.0f\n',length(SpecHead.RangeKm));
fprintf(fid,'%%SpectraDopplerCells: %4.0f\n',SpecHead.PC.SpeclengthPnts);
fprintf(fid,'%%TableType: LLUV RDL_HFR_DP_LERA_v1\n');   %started 3/2017

%%%% original version with all field following COS methods
%fprintf(fid,'%%TableColumns: 18\n');
%fprintf(fid,'%%TableColumnTypes: LOND LATD VELU VELV VFLG ESPC ETMP MAXV MINV ERSC ERTC XDST YDST RNGE BEAR VELO HEAD SPRC\n'); 

%%% cut columns 5 6 7 8 9 from original, as per Otero instructions
%fprintf(fid,'%%TableColumns: 13\n');
%fprintf(fid,'%%TableColumnTypes: LOND LATD VELU VELV ERSC ERTC XDST YDST RNGE BEAR VELO HEAD SPRC\n'); 

%%% added back column 6, the spatial error estimate, to conform to standards
fprintf(fid,'%%TableColumns: 14\n');
fprintf(fid,'%%TableColumnTypes: LOND LATD VELU VELV ESPC ERSC ERTC XDST YDST RNGE BEAR VELO HEAD SPRC\n'); 

%fprintf(fid,'%%TableRows: %5.0f\n',length(find(isnan(Rclean.RadComp(:,PROC_CONST.iwhich_ave))==0)));
fprintf(fid,'%%TableRows: %5.0f\n',length(ikeep));

fprintf(fid,'%%TableStart:\n');
%%% for 18 columns
%fprintf(fid,'%%%%   Longitude   Latitude    U comp   V comp  VectorFlag    Spatial    Temporal     Velocity    Velocity  Spatial  Temporal X Distance  Y Distance   Range   Bearing   Velocity  Direction   Spectra\n');
%fprintf(fid,'%%%%     (deg)       (deg)     (cm/s)   (cm/s)  (GridCode)    Quality     Quality     Maximum     Minimum    Count    Count      (km)        (km)       (km)    (True)    (cm/s)     (True)    RngCell\n');

%%% for 13 columns
% fprintf(fid,'%%%%   Longitude   Latitude    U comp   V comp   Spatial  Temporal X Distance  Y Distance   Range       Bearing   Velocity  Direction   Spectra\n');
% fprintf(fid,'%%%%     (deg)       (deg)     (cm/s)   (cm/s)    Count    Count      (km)        (km)       (km)       (True)    (cm/s)     (True)    RngCell\n');

%%% for 14 columns
fprintf(fid,'%%%%   Longitude   Latitude    U comp   V comp     Spatial Error  Spatial  Temporal X Distance  Y Distance   Range       Bearing   Velocity  Direction   Spectra\n');
fprintf(fid,'%%%%     (deg)       (deg)     (cm/s)   (cm/s)       (cm/s)        Count    Count      (km)        (km)       (km)       (True)    (cm/s)     (True)    RngCell\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% file format
%1 Longitude (deg) 
%2 Latitude   (deg)
%3 U comp (cm/s)   
%4 V comp (cm/s)
% cut %5 VectorFlag   (GridCode)  
%  %6 Spatial Quality   
% cut %7 Temporal   Quality     
% cut %8 Velocity    Maximum   
% cut %9 Velocity    Minimum  
%10  Spatial Count   
%11 Temporal Count     
%12 X Distance  (km)     
%13 Y Distance (km)     
%14    Range  (km)   
%15    Bearing    (True)   
%16    Velocity   (cm/s)    
%17  Direction     (True)   
%18  Spectra RngCell

%%%
for j=1:length(ikeep);
    s=Rclean.RangeBearHead(ikeep(j),:);    
   R=[Rclean.LonLat(ikeep(j),:)'   ;
       Rclean.U(ikeep(j),PROC_CONST.iwhich_ave)   ;
       Rclean.V(ikeep(j),PROC_CONST.iwhich_ave)    ;
      % PROC_CONST.COS_Flag;
       %%%Rclean.Flag(ikeep(j),iwhich_ave)    ;
       Rclean.Error(ikeep(j),PROC_CONST.iwhich_ave)  ;
      % PROC_CONST.Temp_qual    ;
      % PROC_CONST.VelocityMaximum   ; 
      % PROC_CONST.VelocityMimimum   ; 
         
       Rclean.Flag(ikeep(j),2)-floor(Rclean.Flag(ikeep(j),2)./100)*100 ;

       PROC_CONST.TemporalCount    ;
%       s(1)*cosd(true2math(s(3)))   ;   
%       s(1)*sind(true2math(s(3)))   ; 
       s(1)*cosd(s(2))   ;   
       s(1)*sind(s(2))   ; 
   
       Rclean.RangeBearHead(ikeep(j),1)
       math2true(Rclean.RangeBearHead(ikeep(j),2))   ;
       Rclean.RadComp(ikeep(j),PROC_CONST.iwhich_ave)    ;
       math2true(Rclean.RangeBearHead(ikeep(j),3)   ) ;
       round(Rclean.RangeBearHead(ikeep(j),1)./(SpecHead.RangeResKm))] ;
       
    if isnan(R(3))==0  %proceed only if good data exists in this entry
i=find(isnan(R)==1);
R(i)=999;

%%%with length(R)==18
%str='   %3.7f  %3.7f  %3.3f\t%3.3f\t %3.0f\t   %3.3f\t%3.3f\t   %3.3f\t %3.3f\t%2.0f\t%2.0f   %2.4f\t%2.4f  \t%2.4f\t %3.1f\t%2.3f\t   %3.1f  \t%2.0f\n' ; 

%%%with length(R)==13
%%%cut columns 5 6 7 8 9 from original
%str='   %3.7f  %3.7f  %3.3f\t%3.3f\t %2.0f\t%2.0f      %2.4f    %2.4f  \t%2.4f  \t %3.1f\t%2.3f\t   %3.1f  \t%2.0f\n' ; 

%with length(R)==14
%%% add back column 6 (spatial ave) to  length(R)==13 estimate
str='   %3.7f  %3.7f  %3.3f\t%3.3f\t %3.3f\t %2.0f\t%2.0f      %2.4f    %2.4f  \t%2.4f  \t %3.1f\t%2.3f\t   %3.1f  \t%2.0f\n' ; 

fprintf(fid,str,R);
%    pause
    end       
end
%%% done with data TABLE
fprintf(fid,'%%TableEnd:\n');
fprintf(fid,'%%%%\n');
t=datestr(now,30);
fprintf(fid,'%%ProcessedTimeStamp: %s %s %s  %s %s %s \n',t(1:4),t(5:6),t(7:8),t(10:11),t(12:13),t(14:15));

%%
%%% finish with the processing steps used...
for ii=1:length(Rclean.ProcessingSteps)
    s=Rclean.ProcessingSteps{ii};
    i=find(s=='_');
    if  isempty(i)==0 & s(i(end)+1)=='v'
        s1=s(1:i(end)-1);
        s2=s(i(end)+2:end);
        fprintf(fid,'%%ProcessingTool: "%s" %s.0\n', s1,s2 );
    else
    end
end
  %%%% this will output the results in a format like:
%ProcessingTool: "RadialMerger" 10.7.1
%ProcessingTool: "SpectraToRadial" 10.9.0
%ProcessingTool: "RadialSlider" 11.2.0
%ProcessingTool: "RadialArchiver" 11.2.2
%ProcessingTool: "AnalyzeSpectra" 10.7.4

%%% all done, get out.
fprintf(fid,'%%End:\n');
fclose(fid);

