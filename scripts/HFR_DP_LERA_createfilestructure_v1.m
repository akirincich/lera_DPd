function  [site_name]=HFR_DP_LERA_createfilestructure_v1(site_name,base_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HFR_DP_LERA_createfilestructure_v1.m
%
% Script looks for and, if not present, creates the file structure needed to 
%  use and save HF radar data files from LERA-type systems.
%
%
%  Version:
%  -v1  created 
%
%
%  Anthony Kirincich
%  WHOI-PO
%  akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%move to data directory
eval(['cd ' base_dir])

%%
adders={'';'_ts';'_radave';'_pics';'_css';'_config';'_radave_lluv'};

for ii=1:length(adders)
    d=dir(['Site_' site_name adders{ii}]);
    if isempty(d)
        eval(['mkdir Site_' site_name adders{ii}])
    end
end


return
%%
