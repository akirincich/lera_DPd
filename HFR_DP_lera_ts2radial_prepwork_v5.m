%HFR_DP_lera_ts2radial_prepwork_v?
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HFR_DP_lera_ts2radial_prepwork_v?
%
% this script sets up the proper working directories
% and file directories to process the site's spectra
% files and isolates the new files in the css folder that 
% need processing.
%
% version:
%  8/10/2017 started from HFR_DP version
%
%  v3 5/2018 allows the difference in the filenames of the MK2 and MK3 ts files
%
% v4  11/2018
% clean up of v3 and streamline of constants and how the pattern and
% radar_header files are loaded.
%
% v5  1/2019  allowed for new measured pattern to be uploaded, changes to
% load_lera_pattern and radar_pattern.  This uses the data found in 
%  lera_time2spectra_forcal_v1 to get the pattern data
%
% Anthony Kirincich
% WHOI PO
% akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set the directories of files
%processing_file_dir=['/Codar/SeaSonde/Data/Processings/HFR_DP/Site_' site_name];
incoming_ts_file_dir=[base_dir '/Site_' CONST.site_name '_ts'];
incoming_spectra_file_dir=[base_dir '/Site_' CONST.site_name '_css'];
incoming_config_file_dir=[base_dir '/Site_' CONST.site_name '_config'];
outgoing_radialmat_file_dir=[base_dir '/Site_' CONST.site_name ''];
outgoing_radialavemat_file_dir=[base_dir '/Site_' CONST.site_name '_radave'];
outgoing_radialavelluv_file_dir=[base_dir '/Site_' CONST.site_name '_radave_lluv'];
outgoing_pics_file_dir=[base_dir '/Site_' CONST.site_name '_pics'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load header and patt files
eval(['cd ' incoming_config_file_dir])
%load site header
confdir=[incoming_config_file_dir '/'];
config_dir=confdir;   %old 
radar_header %radar_header is used here as well as on the site computer to load this info for the compression step


%%% transfer constants that won't change %%%
CONST.c=299792458;  %m/s speed of light in a vaccum
%     %    c=299705000;   % in air   
%     %    c=299644310;   % in air based on getting v=0 at ibragg (difference of 2 mm/s from first value....)
%number of antenna elements
CONST.N=RC.NANT;

%load site pattern file
%%% use generalized pattern loading program.
%patt=load_lera_pattern_v1(CONST,RC);  %includes ideal for both lpwr and nwtp
patt=load_lera_pattern_v2(CONST,RC);  %includes ideal for both lpwr and nwtp

%now go back to the working directory
eval(['cd ' working_dir])

%%% note that patt.A for the lera should be 1 deg resolution, 360 deg in
%%% span and have the ends of the phi array set to be in the back of the
%%% radar, (i.e. onshore where the signal returns are low).  The last
%%% quality ensures that the peak finder will avoid errors related to peaks
%%% on/near the boundary.  See pattern file and the spectra1radialmetrics file 
%%% for a description of the orientation

%sart HEAD from parts of RC
HEAD=[];
HEAD.SiteName=CONST.site_name;
HEAD.lera_SiteName=RC.SiteName;
HEAD.Site_loc=CONST.Site_loc;
HEAD.Musicparams123=[40 20 2];   %from LPWR
HEAD.ProcessingSteps=[];
HEAD.Bearing=patt.Array_bearing;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%% adjuster for matching pattern
%a=-5;
a=0;
 
%%% The pattern bearings are CCW (counter-clockwise) degrees referenced 
%%% from the antenna bearing. The antenna bearing is found in Header.txt
%%% and is (CW) clockwise degrees from true North. See the File_RadialSetups guide.
%%% for details.
%adjust angles to be in 'true'
patt.angles=patt.Array_bearing+a-patt.angles;
i=find(patt.angles>360); patt.angles(i)=patt.angles(i)-360;
i=find(patt.angles<0); patt.angles(i)=patt.angles(i)+360;
 
%temporarily store full pattern
temp_angles=patt.angles; temp_A=patt.A;

%chop pattern to be only the offshore components wrapping to land,
%use CONST.pattern_ends=[ ] here
i=find(patt.angles>CONST.pattern_ends(1) & patt.angles<CONST.pattern_ends(2));
patt.angles=patt.angles(i);
patt.A=patt.A(:,i);


%%
% make bearing limits for radave files...
% this is harder than it looks to do it a systematic way for all types of patterns
%  WHY? patt.angles starts on one side of the pattern and progresses to the other
%  regardless if there is a discontinuity in the true angles (i.e. crossing 0/360 
%  with jumps).  
% 
% The goal is to make a radave angle array that matches the 
%  direction and sense of rotation of the bearing as patt.angle, but has 
% standardized angles.  This will limit the radial concat issues one might have
% later if the pattern or Antbear is changed.

% %%%%%% old way...this doesn't work if the site wraps around True North
% %%%%%% but is not full coverage
% bearing_limits=[min(patt.angles)+2 max(patt.angles)-2];
% %establish 5deg bearing angle centers and ranges
% BearT=[bearing_limits(1):CONST.bearing_width:bearing_limits(2)];
% BearT_ends=[BearT-CONST.bearing_width/2 BearT(end)+CONST.bearing_width/2];

%%%%%% alternative...define as having centers such that no width spans 0/360
stock_Bear=CONST.bearing_width/2:CONST.bearing_width:360-CONST.bearing_width/2;
db=abs(nanmedian(diff(patt.angles)));
BearT=nan*ones(size(patt.angles));
for ii=1:length(patt.angles)
    [s,i]=sort(abs(stock_Bear-patt.angles(ii))); %find the nearest ave bearing
    if s(1)<= CONST.bearing_width/2 ; %if it is less than CONST.bearing_width/2 away
        BearT(ii)=stock_Bear(i(1));
    end
end

%%% now march through BearT and save only the first value, if there are
%%% more values than CONST.bearing_width/2, so, 3 or more for a 
u=unique(BearT);  
% can't just use unique(BearT,'stable') as also wish to kick out the ones
% that have less than CONST.bearing_width/db/2 values.
for ii=1:length(u)
i=find(BearT==u(ii));
if length(i)>CONST.bearing_width/db/2 % /db is included to account for APMs that don't have increments of 1
    BearT(i(2:end))=nan;
else %cut them all
    BearT(i(1:end))=nan;
end
end
%%% condense
BearT=BearT(~isnan(BearT));
%%% make ends.
BearT_ends=[BearT-CONST.bearing_width/2; BearT+CONST.bearing_width/2];
%%% fix to correct the 0 end if have a Bear? on the high side.
i=find(min(BearT_ends)==0 & [diff(BearT_ends)> CONST.bearing_width] );
if isempty(i)==0 %need to change a 0 to 360
    j=find(BearT_ends(:,i)==0); BearT_ends(j,i)=360;
end
%%% convert to math
Bear=true2math(BearT);
Bear_ends=true2math(BearT_ends);
i=find(min(Bear_ends)==0 & [abs(diff(Bear_ends))> CONST.bearing_width] );
if isempty(i)==0 %need to change a 0 to 360
    j=find(Bear_ends(:,i)==0); Bear_ends(j,i)=360;
end

%%% cut average bearings to be +/-90deg from shore normal
%%% because, while I want the DF to use the full pattern, I don't want the
%%% averaging to use it.
i=find(BearT>=CONST.Site_bounds(1) & BearT<=CONST.Site_bounds(2));

patt.Bear=Bear(i);
patt.Bear_ends=Bear_ends(:,i);
patt.BearT=BearT(i);
patt.BearT_ends=BearT_ends(:,i);

 %%% restore full pattern ... if you wish to output the RM using the full
 %%% pattern.
 % patt.angles=temp_angles; patt.A=temp_A;

%%% make a plot of the pattern, if wanted
%if CONST.goplot(2)==1
    figure(2); clf
%    p=polar(true2math(patt.angles)*pi./180,abs(patt.a1),'r'); hg
    p=polar(true2math(patt.angles)*pi./180,abs(patt.A(1,:)),'r'); hg
     set(p,'linewidth',2)
%     p=polar(true2math(patt.angles)*pi./180,abs(patt.a2),'b');
%     set(p,'linewidth',2)
     p=polar([true2math(HEAD.Bearing)*pi/180]*[1 1],[0 1],'k'); set(p,'linewidth',2);
     title([patt.Site_name ' ' CONST.which_patt ' pattern (converted from DegT to math)']);
    %%% add the azimuths of the radave product, for show
    p=polar(Bear*pi./180,ones(size(Bear)),'kx');
    
    figure(3); clf;
    for iii=1:CONST.N
        subplot(2,4,iii)
        plot(patt.angles,real(patt.A(iii,:)),'b.'); hg; axis([0 360 -1.5 1.5])
        plot(patt.angles,imag(patt.A(iii,:)),'r.'); hg;
%         if strcmp(CONST.which_patt,'meas')==1;
%              plot(patt.angles,real(patt.A_ideal(iii,:)),'g.'); hg; axis([0 360 -1.5 1.5])
%             plot(patt.angles,imag(patt.A_ideal(iii,:)),'m.'); hg;
%         end
     end
    
    
%end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get list of spectra in ts folder
%f=dir([incoming_spectra_file_dir '/CSS*' site_name '*cs*']);
%f=dir([incoming_ts_file_dir '/TS_' RC.SiteName '*.mat']);
f=dir([incoming_ts_file_dir '/*_' RC.SiteName '*.mat']);
fnames={}; fdates=[];
for i=1:length(f)
    fnames(i)={f(i).name};
    aa=char(fnames(i));
%%% for the wera file structure
if strcmp(RC.radar_type,'MK2')==1
    j=find(aa=='_');
    yyyy=str2num(aa(j(2)+1:j(2)+4)); yday=str2num(aa(j(2)+5:j(2)+7)); hh=str2num(aa(j(2)+8:j(2)+9)); mm=str2num(aa(j(2)+10:j(2)+11));
else
    %%% for the lera filename structure
    yyyy=str2num(aa(1:4)); yday=str2num(aa(5:7)); hh=str2num(aa(8:9)); mm=str2num(aa(10:11));
    %year is between 2 and 3, month 3 and 4, day 4 and 5, time 5 and 6
%   fdates(i)=datenum(2000+str2num(aa(j(2)+1:j(2)+2)),str2num(aa(j(3)+1:j(3)+2)),str2num(aa(j(4)+1:j(4)+2)),str2num(aa(j(5)+1:j(5)+2)),str2num(aa(j(5)+3:j(5)+4)),0);
end
fdates(i)=datenum(yyyy,1,0) + yday + hh/24 + mm/(24*60);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get list of RM files in RM  folder
f=dir([outgoing_radialmat_file_dir '/RM*' CONST.site_name '*mat']);
fmnames={}; fmdates=[];
for i=1:length(f)
    fmnames(i)={f(i).name};
    aa=char(fmnames(i));
    j=find(aa=='_');
    %year is between 2 and 3, month 3 and 4, day 4 and 5, time 5 and 6
    fmdates(i)=datenum(str2num(aa(j(2)+1:j(2)+4)),str2num(aa(j(3)+1:j(3)+2)),str2num(aa(j(4)+1:j(4)+2)),str2num(aa(j(5)+1:j(5)+2)),str2num(aa(j(5)+3:j(5)+4)),0);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% identify which of the files in fnames to process

%%% Method 1: choose files based on the file time
if CONST.files_to_process_method==1;
    iwhich=find(fdates>=CONST.files_to_process_dates(1) & fdates<CONST.files_to_process_dates(2) );
    
%%% Method 2: compare css files to RM files, and process new files
elseif CONST.files_to_process_method==2;    
    if isempty(fmdates)==1
        iwhich=1:length(fdates);
    elseif isempty(fmdates)==0
        iwhich=find(fdates>fmdates(end));
    end
end


%%% return the files to process to running script in the fname and fdate
%%% arrays
fnames=fnames(iwhich);
fdates=fdates(iwhich);


