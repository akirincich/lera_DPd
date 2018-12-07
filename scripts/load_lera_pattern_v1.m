function patt=load_lera_pattern_v(CONST,RC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function patt=load_lera_pattern_v(CONST);
%
% This function works to produce a 'beam' pattern structure for a 
%  'lera' type HFR system with 8 antenna channels.  The antenna
%  can be arbitrarily spaced and arranged, but the map that describes the
%  spacing and orientation of the array must be described here in terms of
%  a matrix of antenna response patterns for each entry of a steering
%  vector (i.e. the bearing of a plane wave impinging on the array.  
%
%  To be clear, the patt structure MUST look like this:
%
% patt.Site_name        %the site name 
% patt.date             %the date of the pattern or the current date if ideal
% patt.Pattern_type     %the 'ideal' or 'meas' flag
% patt.mapper           % an CONST.N x 3 matrix with a map from the antenna to the rx channel
% patt.mapper_what      % usually '(1) ideal_pattern_ant# (2) array channel  (3) map to transform ideal # to channel 1-8';
% 
% patt.angles=psi;      %reports the angle the target is coming from, not going toward.
%                       %  This is the response angle, not the steering angle, 
%                       %    and is defined here as the direction a wave-
%                       %   form would be coming from as measured by the array
%                       %   design, following van Trees.
% 
% patt.A                %an CONST.N x patt.angle matrix of the antenna pattern response relative to
%                       % a center antenna of the array, and relative to the coo,  thus need to know how much to 
%                       %rotate this by to put into earth coordinates
% 
% patt.Array_bearing    % as defined, the direction the pattern x-axis is pointing
% 
% patt.Array_bearing_units  %'Degrees cw from True North';
% patt.sense            %a descriptor of the pattern as defined in radar_pattern.m
% 
% 
%   As a example:
%   using the RX pattern response when arranged as a rectangular
%   grid.  To create the 'ideal' pattern we define:
%
%   N     - number of elements in x direction, 
%   M     - number of elements in y direction 
%   phi   - steering angle in math coordinates (in degrees relative to x (eg ccwE))*
%   counterclockwise of east (if the x axis of the instrument array was 
%    pointed east, which it never is)
%   
%  NOTE: the 'steering angle' as defined here following van Trees is
%  the direction a waveform is going TOWARDS, not the direction it is
%  COMING from.   In HFR processing, we use the array as a RX array and are 
%   solely interested in the direction a target waveform is coming from 
%  and how that would map onto the antenna array.
%
%   THUS, to convert this into a useful map for an antenna pattern
%     measurement, and to be able to compare the ideal to measured ant
%     response patterns we first follow van Trees to develop the antenna
%     beam response pattern, and then define the 'response angle' as
%     psi=phi+180.  The response angle is returned within patt as psi to define
%     the antenna response in a format consistent with our HFR needs.
%
%  following van Trees in the context of our rectangular array:
%
% ^
% |          x   o   o
% y          o   o   o
%            o   o   o
%
%                x ->
%
% As defined here:
%     (1) the 'bearing' of the array should be the direction (heading) of the x axis.
%     (2) the 'response angle' is defined as ccw of the positive x-axis through 
%            the array
%
%
% Version
%   v1  created from lera_pattern_work_v4.m and adjusted for loading ideal
%   or measured pattern generalities
%
% by
% Anthony Kirincich
% WHOI PO
% akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(CONST.which_patt,'Ideal')==1

%%%%%%% set up outgoing variable
patt=[];
patt.Site_name=CONST.site_name;
%%% pull out the date of the pattern calculation.
patt.date=date;
patt.Pattern_type=CONST.which_patt;
    
%%% since each radar can have the rx antennas orientated in slightly
%%% different configurations, call the radar_pattern.m file to fill in the
%%% radar specific variables that would fill out the patt structure.
%%%
%%% In general, we need to get:
%%%      psi
%%%      A 
%%%      mapper 
%%%

radar_pattern;  %%%% get site specific information here

%so the cable number that matches the ant_num is
patt.mapper=mapper;
patt.mapper_what=mapper_what;
patt.A=A(mapper(:,3),:);

%%% package this up and add the ancillary data needed.
patt.angles=psi0;   %report the angle the target coming from, not going toward.

patt.Array_bearing=Array_bearing;
patt.Array_bearing_units='Degrees cw from True North';
patt.sense=sense;


%%

elseif strcmp(CONST.which_patt,'Measured')==1
    %if there is a MeasPatt file in the *_config directory, find it and
    %load it. 
    %if there are multiple MeasPatt files in the directory, load the most
    %recent one.
    
    %look for MeasPatt files in the current directory, should be the config dir
    d=dir('Meas_Pattern.mat');
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

