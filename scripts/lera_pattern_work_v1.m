function patt=lera_pattern_work_v(site_name,which_patt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function patt=lera_pattern_work_v(site_name,which_patt);
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear 
%cd ~/Matlab/working/LERA


% N     - number of elements in x direction, 
% M     - number of elements in y direction 
% phi   - steering angles (in degrees relative to x (eg ccwE))*
%   i.e. math coordinates
%   counterclockwise of east (if the x axis (of the instrument pointed east
%
%
% but here this is ccw of the positive x-axis through the array
% ^
% |          o   o
% y          o   o   o
%            o   o   o
%
%                x ->
%

if strcmp(which_patt,'ideal')==1
    
%phi=180:1:360;
phi=-179:1:180;
%phi=0:1:359;
th=90;  %the vertical orientation of the array slice, 90 deg is horizontal
M=3; N=3;

% Define dx/lambda and dy/lambda, that is, normalized by lambda
%dx = 0.5; dy = 0.5; 

%actual spacing of array is .5*lambda on diagonal, so 
%(.5*l)^2 = 2 * (x*l)^2
%x*l = sqrt((.5*l)^2 /2)
%x = sqrt( (.5*l)^2 /2)/l;   
%  x =  sqrt(  (.5^2 * l^2) / (l^2 * 2) = sqrt ( .5^2/2 )
%x = 0.3535  

dx = 0.3535; dy = 0.3535; 


% See Van Trees, 2002, pg 235, 240, 249

% Make grid of M and N, indexed from zero, M columns
% % this seems to put the phase center at a corner
% [m,n] = meshgrid(0:M-1,0:N-1);
%
% Try to put the phase center in the middle
[m,n] = meshgrid(  (0:M-1) - mean(0:M-1)  ,   (0:N-1) - mean(0:N-1)   );
%n=flipud(n);
%
% m = -1     0     1
%     -1     0     1
%    -1     0     1
%n = -1    -1    -1   or after flipud:    1     1     1
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
A = exp( 1i.* ( (n .* psi_x)  + (m .* psi_y) ) );
%
%cut A down to 8 channels as the top right corner of the array is the missing antenna
%[ m(:,1) n(:,1) ] =     -1    -1      ant_num   1
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
%    x,y    (1)  -1,-1    (4)   0,-1   (7)   1,-1
%           (2)  -1,0    (5)   0,0   (8)   1,0
%           (3)  -1,1   (6)   0,1  (9)   1,1

sense= {'Sense of rotation is CW is negative, CCW is positive, i.e. math';
 ' orientation is about the upwards y-axis or towards ant 4,';
 ' i.e.  a hit at at phi=0 means the scatterer is comming from the direction of';
 '  ant4, (as defined, phi is the direction along the crests so we rotate 90 deg to get the';
 '  wavenumber, as done here)'};

% so point 9 is the one to cut.
%A=A(1:8,:);  %wrong?  7/29

% NO! point 7 is the one to cut
%A=A([1:6 8:9],:);

% NO! NO! point 1 is the one to cut, reorient the whole array to 
A=A([2:9],:);

%%% package this up and add the ancillary data needed.
patt=[];
patt.Site_name=site_name;
%patt.A=A;
%A is the antenna pattern response relative to the center line of the
%antenna array, thus need to know how much to rotate this by to put into
%earth coordinates
patt.Array_bearing=[147]; % based on bearing of 5 to 3 ant line. (see below)

%the compass bearing (CW of true north) of the positive x direction of the 
% antenna pattern relative to the center of the array.

%%% pull out the date of the pattern calculation.
patt.date=date;
patt.Pattern_type='ideal';
patt.angles=phi(1,:);
 %Need lambda for this?
patt.fc=16.15;  % center frequency (in MHz)

%%%%%%%%%%
%%% a look at the Rx back end  %%%   
%         channel number  8  7  6  5  4  3  2  1
% 6/1/2017  cable number  9  8  7  6  5  4  3  2    with cable 1 as a spare
% 7/18/2017 cable number  9  8  7  6  5  4  1  2

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

%  so: ideal actual
mapper=[ 1  1; %8;
         2  4; %3;
         3  5; %4;
         4  2; %7;
         5  3; %2;
         6  6; %1;
         7  7; %6;
         8  8]; %5];
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
% mapper=[ 1  2;
%          2  4;
%          3  5;
%          4  1;
%          5  3;
%          6  6;
%          7  7;
%          8  8];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      

 [s,i]=sort(mapper(:,2),'ascend'); mapper(:,3)=i;

%so the cable number that matches the ant_num is
patt.mapper=mapper;
patt.mapper_what='(1) ideal_pattern_ant# (2) array channel  (3) map to transform ideal # to channel 1-8';

patt.A=A(mapper(:,3),:);
patt.sense=sense;
%%%%% working parts %%%%

elseif strcmp(which_patt,'meas')==1
    
    asdfasdfasdf
    
end %strcmp

%%




% 
% %cut A down to 8 channels as the top right corne of the array is the missing antenna
% %[ m(:,1) n(:,1) ] =     -1    -1      ant_num   1
% %                        -1     0                2
% %                        -1     1                3
% %                         0    -1                4
% %                         0     0                5
% %                         0     1                6
% %                         1    -1                7
% %                         1     0                8
% %                         1     1                9
% %
% % or:
% %    x,y    (1)  -1,-1    (4)   0,-1   (7)   1,-1
% %           (2)  -1,0    (5)   0,0   (8)   1,0
% %           (3)  -1,1   (6)   0,1  (9)   1,1
% 
% sense={'Sense of rotation is CW is negative, CCW is positive';
% ' orientation is about the upwards y-axis or towards ant 4,';
% ' i.e.  a hit at at phi=0 means the wave is comming from the direction of';
% '  ant4, phi is the direction along the crests so rotate 90 deg to get the';
% '  wavenumber'};
% 
% % so point 9 is the one to cut.
% %A=A(1:8,:);  %wrong?  7/29
% 
% % NO! point 7 is the one to cut
% %A=A([1:6 8:9],:);
% 
% % NO! point 1 is the one to cut, reorient the whole array to 
% A=A([2:9],:);
% 
% %%% package this up and add the ancillary data needed.
% patt=[];
% %patt.A=A;
% %A is the antenna pattern response relative to the center line of the
% %antenna array, thus need to know how much to rotate this by to put into
% %earth coordinates
% patt.Array_bearing=[90]; 
% 
% %the compass bearing (CW of true north) of the positive x direction of the 
% % antenna pattern relative to the center of the array.
% 
% %%% pull out the date of the pattern calculation.
% patt.date=date;
% patt.Pattern_type='ideal';
% 
%  %Need lambda for this.
% patt.fc=16.25;  % center frequency (in MHz)
%  
% % the  nwt cable grid is:
% %  5  2
% %  4  3  6
% %  9  8  7
% % 
% % which maps into channels
% %  4  1       or really    5  6
% %  3  2  5              1  2  7
% %  8  7  6              4  3  8
% %  so: ideal  actual
% % mapper=[ 1  4; %8;
% %          2  3; %3;
% %          3  8; %4;
% %          4  1; %7;
% %          5  2; %2;
% %          6  7; %1;
% %          7  5; %6;
% %          8  6]; %5];
% mapper=[ 1  1; %8;
%          2  4; %3;
%          3  5; %4;
%          4  2; %7;
%          5  3; %2;
%          6  6; %1;
%          7  7; %6;
%          8  8]; %5];
% 
%  [s,i]=sort(mapper(:,2),'ascend'); mapper(:,3)=i;
% 
% %so the cable number that matches the ant_num is
% patt.mapper=mapper;
% patt.mapper_what='(1) ideal_pattern_ant# (2) array channel  (3) map to transform ideal # to channel 1-8';
% 
% patt.A=A(mapper(:,3),:);
% patt.sense=sense;
% %%%%% working parts %%%%
% %%


% what does this do?
%B = beam_pattern(A, find(phi==0) );


%%
%  %okay, now, what if I have the exact dimensions of the array?  How would I
%  %get a version of A?
%  
%  %Need lambda for this.
%  lambda=16.25;  % in MHZ
%  
% % the  nwt cable grid is:
% %  5  2
% %  4  3  6
% %  9  8  7
% %so the cable number that matches the ant_num is
% cable_num= [ 9 4 5 8 3 2 7 6];
% channel_num=cable_num-1;
% 
% cxy=[ 9  0     0;
%       4  0    21.3;
%       5  0    21.3+21.8;
%       8 3 2 7 6];
% 
% 
% %  start stop meas
% ssm=[  ]
% 
%  then divide m,n by lambda.
%  
%  Need pointing angle of array as well.
%  
%  not too hard... right?
%      
% 
% 
% then take A into direction finding, along with cospectra matrix...
% 
 
 
 
 
 
 
 
