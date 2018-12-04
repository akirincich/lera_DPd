function patt=lera_pattern_work_v(site_name,CONST);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function patt=lera_pattern_work_v(site_name,which_patt);
%
% This function works to produce a 'beam' pattern structure for a 
%  'lera' type HFR system with 8 antenna channels arranged as a rectangular
%  grid.  To create the 'ideal' pattern we define
%
% N     - number of elements in x direction, 
% M     - number of elements in y direction 
% phi   - steering angle in math coordinates (in degrees relative to x (eg ccwE))*
%   counterclockwise of east (if the x axis of the instrument array was 
%    pointed east, which it never is)
%   
%  NOTE: the 'steering angle' as defined here following van Trees is
%  the direction a waveform is going TOWARDS, not the direction it is
%  COMING from.   In HFR processing, we use the array as a RX array and are 
%   solely interested in the direction a target waveform is coming from 
%  and how that would map onto the antenna array.
%     THUS, 
%     to convert this into a useful map for an antenna pattern
%     measurement, and to be able to compare the ideal to measured ant
%     response patterns we first follow van Trees to develop the antenna
%     beam response pattern, and then define the 'response angle' as
%     psi=phi+180.  The response angle is returned within patt as phi to define
%     the antenna response in a format consistent with our HFR needs.
%
%  following van Trees in the context of our array:
%
% but here this is ccw of the positive x-axis through the array
% ^
% |          x   o   o
% y          o   o   o
%            o   o   o
%
%                x ->
%
% Thus the 'bearing' of the array should be the direction of the x axis.
%  (a change from vers. 1 of this script.
%
%   see below for additional details
%
% . Version
% . v1: created, wrong, and edited 7/2017 and 8/2017
%
% .  v2:  9/2017
% .    adjusted to current for of array design 
% .        corrected the n,m synatacs from van Trees
% .    added the need for the CONST for advanced plotting of examples (and
%        to ensure non overlap with v1, (b/c orientation, bearing angles
%        are changed.
%     corrected ideal set up to focus end result on the direct a waveform 
%      was coming from, not going toward.
%   
%    reporting patt.angle as psi not phi.
%
%    NWT array was reset in July 17,2017 because of cable break on the
%    original cable 3, this script accounts for this change
%
%   9/18/2017 . also added the ability to load a measured pattern file.
%               measured patterns for lera systems are formed by starting
%               with the ideal pattern result and comparing this to the
%               combined gps and timeseries data.
%             see HFR_DP_lera_calwork_v08112017.m for how this was done for
%             NWPT.   other sites will likely be different.
%
% v3 3/2018   uses a different way to do phi,psi formulation
%         %%%  not any different from v2 and has some mistakes in the phi, psi formulation.
%          %%% DO NOT USE %%%
%
%   v4 5/2018   fixed for multiple sites, differentiated by the placement
%                  of the antenna
%
%
% by
% Anthony Kirincich
% WHOI PO
% akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%%% for testing
% clear 
% cd ~/Matlab/working/LERA
% 
% site_name='NWPT';
% which_patt='ideal';
%%%%%%%%%%%%%%


if strcmp(CONST.which_patt,'ideal')==1

%%%%%% set array variables %%%%
%phi=180:1:360; 
phi=1:1:360;  % the steering angle, defined following van Trees.
psi= phi-180;   %  the response angle, defined here as the direction a wave-
                  %   form would be coming from as, measured by the array
                  %   design, following van Trees.
%%% both phi and psi are in math coordinates!!! %%%
 
th=90;  %the vertical orientation of the array slice, 90 deg is horizontal
M=3; N=3;
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%  define element spacing.
% Define dx/lambda and dy/lambda, that is, normalized by lambda
dx = 0.3535; dy = 0.3535;   %for lera system, ideal
% See Van Trees, 2002, pg 235, 240, 249

% Make grid of M and N, indexed from zero, M columns
% % this seems to put the phase center at a corner
% [m,n] = meshgrid(0:M-1,0:N-1);
%
% Try to put the phase center in the middle
%[m,n] = meshgrid(  (0:M-1) - mean(0:M-1)  ,   (0:N-1) - mean(0:N-1)   );
[n,m] = meshgrid(  (0:N-1) - mean(0:N-1) , (0:M-1) - mean(0:M-1)  );   %switch m,n to match van trees

m=flipud(m);  %flipup needed to get the same sense as van trees fig 4.7
%
% n = -1     0     1
%     -1     0     1
%    -1     0     1
% m = -1    -1    -1   or after flipud:    1     1     1
%     0     0     0                       0     0     0       
%     1     1     1                      -1    -1    -1 
%     

% vectorize these, down each column, row 1 to m
m = m(:);
n = n(:);

% now expand to size of phi
m = repmat(m,1,length(phi));
n = repmat(n,1,length(phi));

% make matrix out of phi
phi = repmat(phi(:)',size(m,1),1);

% eqn 4.2, 4.3 with dx defined in wavelengths
psi_x = 2.*pi.*dx.*sind(th).*cosd(phi);
psi_y = 2.*pi.*dy.*sind(th).*sind(phi);


% Compute the matrix of array manifold vectors (eqn 4.50, 4.53)
A = exp( 1i.* ( (n .* psi_x)  + (m .* psi_y) ) );  %follows van trees, working from equ 4.1 to 4.50

% for the arbitrary array orientation, shown above...
%[ n(:,1) m(:,1) ] =     -1    -1      ant_num   1
%                        -1     0                2
%                        -1    1                3
%                         0     -1                4
%                         0     0                5
%                         0     1               6
%                         1     -1                7
%                         1     0                8
%                         1    1                9
%
% or:
%    x,y    (1)  -1,1    (4)   0,1   (7)   1,1
%           (2)  -1,0    (5)   0,0   (8)   1,0
%           (3)  -1,-1   (6)   0,-1  (9)   1,-1


sense= {'Sense of rotation is CW is negative, CCW is positive, i.e. math';
 ' orientation is about the rightward x-axis or towards ant 8,';
 ' i.e.  a hit at phi=0 (psi=180) means the scatterer (and the wavenumber vector) is going from ant 2 towards ant 8';
 '  (as defined, the crest are normal to the wavenumber vector)';
  ' THUS';
  'for a scattered waveform coming from offshore of ant 9 and traveling toward ant 1:';
  'phi= 135 and psi = -45 ';
  ' i.e. psi is still in math coordinates about the bearing direction of the array (the pos x axis)'};

%%
if CONST.goplot(2)==1;
    %%%%% plot some examples for reference
    tt=[360 45 90 135];
    figure(10); clf;
    [np,mp] = meshgrid(  (0:N-1) - mean(0:N-1) , (0:M-1) - mean(0:M-1)  );   %switch m,n to match van trees
    mp=flipud(mp);  %flipup needed to get the same sense as van trees fig 4.7
    
    for ii=1:length(tt)
        i=find(phi(1,:)==tt(ii));
        x=[real(A(1:3,i))'; real(A(4:6,i))'; real(A(7:9,i))']';
        y=[imag(A(1:3,i))'; imag(A(4:6,i))'; imag(A(7:9,i))']';
        subplot(1,length(tt),ii);
        contourf(np,mp,x,'linecolor','none'); caxis([-1 1]); hold on; grid on;
        contour(np,mp,y,[0 0],'linecolor','w','linestyle','-');
        contour(np,mp,y,[.1:.1:1],'linecolor','k','linestyle','-');
        contour(np,mp,y,[-1:.1:-.1],'linecolor','k','linestyle','--');
        title(['incoming waveform for phi=' num2str(phi(1,i))])
x
y
%pause 
    end
    disp('note that contour is messing up the offangle parts because of interpolation')
end

%%

%%% as described here, in the lera array with the orientation given,
%%% point 1 is the one to cut, reorient the whole array to 
A=A([2:9],:);

%%% package this up and add the ancillary data needed.
patt=[];
patt.Site_name=site_name;
%%% pull out the date of the pattern calculation.
patt.date=date;
patt.Pattern_type='ideal';

%patt.angles=phi(1,:);
patt.angles=psi;   %report the angle the target coming from, not going toward.

 %Need lambda for this?
%patt.fc=16.15;  % center frequency (in MHz)

%patt.A=A;
%A is the antenna pattern response relative to the center line of the
%antenna array, thus need to know how much to rotate this by to put into
%earth coordinates
patt.Array_bearing_units='Degrees cw from True North';
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch site_name

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'NWTP'    %nantucket  LERA, installed on 6/2017

%%% update for new orientation of math here, 
patt.Array_bearing=[147+90]; % based on bearing of 5 to 3 ant line. (see below)
% .                            and adjusted to line with x axis of ant
%                              array, which should be along ants 1 to 7

%the compass bearing (CW of true north) of the positive x direction of the 
% antenna pattern relative to the center of the array.

%%%%%%%%%%
%%% .     a look at the Rx back end  %%%   
%         channel number  8  7  6  5  4  3  2  1
% 6/1/2017  cable number  9  8  7  6  5  4  3  2    with cable 1 as a spare
% 7/18/2017 cable number  9  8  7  6  5  4  1  2    with no spare as cable 3 was cut by a mower

%%%%%% from  6/1 to 7/17/2017 %%%%%%%
% the  nwt cable grid is:
%  5  2     or flipped    6   7
%  4  3  6             2  3   8
%  9  8  7             5  4   9
% 
% which maps into channels
%  4  1       or flipped   5  6
%  3  2  5              1  2  7
%  8  7  6              4  3  8

% %  so: ideal actual
% mapper=[ 1  1; %8;
%          2  4; %3;
%          3  5; %4;
%          4  2; %7;
%          5  3; %2;
%          6  6; %1;
%          7  7; %6;
%          8  8]; %5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% from 7/18/2017 to present %%%%%%
% the  nwt cable grid is:
%  5  2     or flipped    6   7
%  4  3  6             1  2   8
%  9  8  7             5  4   9
% 
% which maps into channels
%  4  1       or flipped   5  6
%  3  2  5              2  1  7
%  8  7  6              4  3  8

% %  so: ideal actual
mapper=[ 1  2;
         2  4;
         3  5;
         4  1;
         5  3;
         6  6;
         7  7;
         8  8];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 [s,i]=sort(mapper(:,2),'ascend'); mapper(:,3)=i;

%so the cable number that matches the ant_num is
patt.mapper=mapper;
patt.mapper_what='(1) ideal_pattern_ant# (2) array channel  (3) map to transform ideal # to channel 1-8';

patt.A=A(mapper(:,3),:);
patt.sense=sense;

disp('NWTP:  Note that the cable mapping changed on 7/18/2017.  Are you using the correct map?')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'LPWR'
        
        %%% update for new orientation of math here, 
patt.Array_bearing=[180+90]; % based on bearing of 5 to 3 ant line. (see below)
% .                            and adjusted to line with x axis of ant
%                              array, which should be along ants 1 to 7

%the compass bearing (CW of true north) of the positive x direction of the 
% antenna pattern relative to the center of the array.

%%%%%%%%%%
%%% .     a look at the Rx back underside  %%%   
%         channel number  8  7  6  5  4  3  2  1
% 5/2018  cable number  8  7  6  5  4  3  2  1   with cable 9 as a spare

%%%%%% from  4/2018 to present  %%%%%%%
% % the  mvy cable grid is:   (from ian email of 5/10/2018
% ~~~~~~~~~ (ocean)
% 9       6       3
% 8       5       2
% 7       4       1
% 
% We use 1 to 8 for channel.
% Dta file should be 
% 1122334455667788
% IQIQIQIQIQIQIQIQIQIQ, my initials

%thus ant and cable pattern is
%     ^   y
%     |  
%     6  3   
%  8  5  2      ---->  x
%  7  4  1       

% %  so: ideal actual
mapper=[ 1  8;
         2  7;
         3  6;
         4  5;
         5  4;
         6  3;
         7  2;
         8  1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 [s,i]=sort(mapper(:,2),'ascend'); mapper(:,3)=i;

%so the cable number that matches the ant_num is
patt.mapper=mapper;
patt.mapper_what='(1) ideal_pattern_ant# (2) array channel  (3) map to transform ideal # to channel 1-8';

patt.A=A(mapper(:,3),:);
patt.sense=sense;

disp('LPWR:  Note that the antenna and cable map is the same, but map is flipped from NWTP')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'HBSR'
        
        %%% update for new orientation of math here, 
patt.Array_bearing=[180+90]; % based on bearing of 5 to 3 ant line. (see below)
% .                            and adjusted to line with x axis of ant
%                              array, which should be along ants 1 to 7

%the compass bearing (CW of true north) of the positive x direction of the 
% antenna pattern relative to the center of the array.

%%%%%%%%%%
%%% .     a look at the Rx back underside  %%%   
%         channel number  8  7  6  5  4  3  2  1
% 5/2018  cable number  8  7  6  5  4  3  2  1   with cable 9 as a spare

%%%%%% from  4/2018 to present  %%%%%%%
% % the  mvy cable grid is:   (from ian email of 5/10/2018
% ~~~~~~~~~ (ocean)
% 9       6       3
% 8       5       2
% 7       4       1
% 
% We use 1 to 8 for channel.
% Dta file should be 
% 1122334455667788
% IQIQIQIQIQIQIQIQIQIQ, my initials

%thus ant and cable pattern is
%     ^   y
%     |  
%     6  3   
%  8  5  2      ---->  x
%  7  4  1       

% %  so: ideal actual
mapper=[ 1  8;
         2  7;
         3  6;
         4  5;
         5  4;
         6  3;
         7  2;
         8  1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 [s,i]=sort(mapper(:,2),'ascend'); mapper(:,3)=i;

%so the cable number that matches the ant_num is
patt.mapper=mapper;
patt.mapper_what='(1) ideal_pattern_ant# (2) array channel  (3) map to transform ideal # to channel 1-8';

patt.A=A(mapper(:,3),:);
patt.sense=sense;

%disp('LPWR:  Note that the antenna and cable map is the same, but map is flipped from NWTP')
disp('LPWR:  Note array is not yet set. ')

        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


elseif strcmp(CONST.which_patt,'meas')==1
    %if there is a MeasPatt file in the *_config directory, find it and
    %load it. 
    %if there are multiple MeasPatt files in the directory, load the most
    %recent one.
    
    %look for MeasPatt files in the current directory, should be the config dir
    d=dir('MeasPatt_*.mat');
    %error if there is no file found.
if isempty(d)==1;
    disp(['Error!  No Measured Pattern file was found in the directory: ' pwd] )
asdfasdf
else
    patt_fns={};    patt_dates=[];
    for ii=1:length(d)
      patt_fns{ii}=d(ii).name  ;
      patt_dates(ii)=d(ii).datenum;
    end   
    %find which of the found files is the newest
    [s,i]=sort(patt_dates,'descend');
    patt_fn=char(patt_fns{i(1)})
    
    load(patt_fn)
end
end %strcmp

