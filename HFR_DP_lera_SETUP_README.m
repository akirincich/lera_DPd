%HFR_DP_SETUP_README.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HFR_DP_SETUP_README.m
%
%  This package contains the master programs to process data from MK-II 
%  or MK-III LERA-type HF Radar systems. 
%  
%  The package requires data from the decimated timeseries produced by: 
%
%        dtacq2mat_compress.m   or equivalent
%
%  on the site computer that collects the timeseries data.  Both the site
%  computer and this offline processing must be coordinated with the same
%  radar setup files:
%
%       radar_header.m
%       dds_prog.out.DD-MON-YYYY.mat   or  radcelf_prog_DD-MON-YYYY.mat
%
%  that describe the sampling characteristics of the radar.
%
%  Within the package, the m-file:
%
%     HFR_DP_master_LERA_XXXX.m     
%
%  where the XXXX denotes the 4-letter radar site code serves as the master
%  script controlling the setup, loading, and processing of the timeseries 
%  data to range and doppler spectra estimates, radial metrics, and radial 
%  average output.
%
%   After reading this readme file, start with the master script to setup 
%  and control your (re)processing.
%
%  This package distributed as is, with no additional warranties, as an open
%  source community resource to enable HFR methodological development, 
%
%  Similar to HFR_Progs, no part of this distrubution can be used for 
%  commerical gain.

%  A series of test master files are included in this release for radar
%  sites LPWR and NWTP.  Raw data for these sites can be found at:
%
%  https://github.com/akirincich/lera_DP_testdata.git
%
%
% v1        December 2018
%
% Anthony Kirincich
% Woods Hole Oceanographic Institution
% akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This README file describes:
%  (1) what HFR_DP is, 
%  (2) what is needed to run HFR_DP, and
%  (3) what HFR_DP is not.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% What HFR_DP_LERA does                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  The HFR_DP package is intended to be used as a platform to reprocess 
%  HFR radar data for direct scientific research and/or advancing HFR meth-
%  ological research.  The program and its structure are based on the 
%  HFR_Progs matlab toolbox (by Kaplan and Cook) and the package contains 
%  versions of key HFR_Progs scripts that have been adapted to provide 
%  additional output or functionality. The steps, outputs, and file structure
%  are taken from the HFR_DP program for COS type radar systems previously 
%  released (www.github.com/akirincich/HFR_DP).
%
%  Testing and validation of this package and its methods are given in:
% 
%  Kirincich, Emery, Washburn, and Flament, "Surface Current Mapping Using 
%  a Hybrid Direction Finding Approach for Flexible Antenna Arrays" 
%  JOAT, (submitted)
%
%  This package uses direction finding methods to determine the emitter 
%   or source locations.
%
%  This package allows for a variety of processing choices to be used and 
%  variable antenna arrays both in terms of number of antennas and in their 
%  arrangement.  The program uses or allows:
%
%  -Image-based First Order Limits (Kirincich, 2017)
%  -MLE, WSF, or MUSIC-type direction finding algorithms
%  -MLD or AIC statistical emitter determination methods
%  -MUSIC_parameter or MUSIC_highest empirical emitter determination methods
%
%  - Radial metrics output (Kirincich et al 2012; de Paola et al 2015;
%     Haynes et al 2017)
%
%  HFR_DP_LERA performs NO temporal averaging of the Radial estimates.  Thus, 
%  the end result is equivalent to COS-processed 'Radial Short' files.
%
%  There numerous types of output products available.  All are available for 
%  each individual spectral estimate submitted for processsing.  Outputs
%  include:
%
%  - Radial Metrics output   (as a Matlab structure within a *.mat file)
%
%  - Radial Average output   (the spatial averaged result of the Radial 
%                             Metrics files, as a matlab structure within
%                             a *.mat file) 
%
%  - LLUV output             (the spatial averaged result of the Radial 
%                             Metrics files, saved in an ASCII text file 
%                             equivalent to COS lluv format and accepted by
%                             NOAA-IOOS HF program for ingestion into the
%                             national network.
%
%  - *.jpg figures of the Spectral results, Radial metrics, and Radial 
%                             averaged results are also available.
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% What HFR_DP Needs                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  HFR_DP requires the following programs to be loaded and operational 
%  on the target machine:
%
%  (1) A MAC or LINUX machine.
%
%  (2) MATLAB (tested on post-2013 releases only)
%
%  (3) The HFR_Progs and M_map packages:
%       find the HFR_Progs at       https://github.com/rowg/hfrprogs.git
%       and the m_map package at    https://www.eoas.ubc.ca/~rich/map.html
%       these packages must be installed and in your path
%
%  (4) The Matlab Image Processing Toolbox is required to run ImageFOLs. 
%
%  (5) A subdirectory file structure for each site_name=XXXX is required, and 
%      must be formulated as such (where ~ denotes some base directory):
%
%      ~/SITE_XXXX_ts          %where the timeseries files will be found
%      ~/SITE_XXXX_css          %where the spectral files will be saved
%      ~/SITE_XXXX_config       %where the header and pattern data will be found
%      ~/SITE_XXXX              %where the resulting 'Radial Metrics' output will go
%      ~/SITE_XXXX_radave       %where the resulting 'Radial Short' will go
%      ~/SITE_XXXX_radave_lluv  %where the resulting ascii file of the 'Radial Short'
%                               %  will go...suitable for transmission to the national archive
%      ~/SITE_XXXX_pics         %where all saved jpg images will go 
%
%      A function: HFR_DP_LERA_createfilestructure_v1.m embedded in the 
%      HFR_DP_master_LERA_XXXX.m file will create those folders which you've 
%      not yet setup.  However, to run master script to process data, at
%      minium, the SITE_XXXX_ts and SITE_XXXX_config folders should be present 
%      and contain data.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% What HFR_DP_LERA doesn't do                      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         HFR_DP does not:
%
%    (1)  Provide significant pre-processing noise cancelation.  
%         Only 2 methods are used to decrease the potential for noise to 
%         influence the radial calculation. 
%         (a) The timeseries2spectra step includes a basic 'hitch' removal feature
%             for a narrow-banded noise sources. 
%         (b) The imageFOL-based delineation of what data is processed has 
%             minimal steps to ignore a spectral area if given characteristics
%             of the spectra are not met, usually due to noise or interference
%         No formal de-striping, etc. is done here.  
%
%    (2)  Look for or account for ships, etc. within the spectral data, 
%
%    (3)  Perform any temporal averaging of the resulting radial velocities.
%         Thus, all resulting data, apart from the radial metrics themselves
%         is equivalent to COS's 'radial shorts' files.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% How to run HFR_DP_LERA                           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    (1) Add the HFR_Progs, M_map, and HFR_DP_LERA packages to your
%        matlab path.
%
%    (2) Create a working directory structure as shown above for each site
%        you intend to reprocess, or download the test data sets also
%        distributed by Kirincich.
%
%    (3) Open, Edit, and Resave:        HFR_DR_master_LERA_XXXX.m 
%        with the information specific to your SITE
%         (i.e. LPWR), directories, and processing choices.
%   
%    (4) Iterate as you need.
%
%    (5) a series of test master files are included in this release for radar
%        sites LPWR and NWTP.  Raw data for these sites can be found at:
%
%        https://github.com/akirincich/lera_DP_testdata.git
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% END OF README
