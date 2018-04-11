function lera_spectra_plot_v1(base_dir,filein);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function lera_spectra_plot_v1(base_dir,filein);
% 
%   A script for plotting of the data from a lera spectral file
%
%   saves jpg to the data/site folder
%
% Anthony Kirincich
%  WHOI_PO
%  akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%% set directories %%%%%%%%%%%%%
%base_dir='~/Matlab/working/LERA/'

%data_dir
data_dir=[base_dir 'lera_data/'];
%script dir
script_dir=[base_dir 'lera_DP/'];
%config_dir
config_dir=[base_dir 'lera_config/'];

addpath(script_dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%% set the file coming in %%%%
% %filein='CSE_nwt_20171701130.mat';
% filein='CSE_nwt_20171712030.mat';

%%% switch back to data dir
eval(['cd ' data_dir]);
%%% load the file
load(filein)

%%

%set up plot of power for all ants
figure(1); clf;
r=SpecHead.RangeKm'*ones(1,length(SpecHead.doppler_velc));
v=ones(length(SpecHead.RangeKm),1)*SpecHead.doppler_velc;

for ii=1:RC.NANT
    eval(['d=data.a' num2str(ii) num2str(ii) ';'])
%    eval(['d=abs(data.a1' num2str(ii) ');'])
    subplot(4,2,ii);

        pcolor(v,r,20*log10(d)); shading flat; hg; caxis([-50 80]); 
        plot([1; 1]*SpecHead.doppler_velc(SpecHead.iFBraggc),[1 1; max(r(:)) max(r(:))],'k')
title(['Antenna ' num2str(ii) ' Ens. Averaged Autospectra'])
if ii==7;  ylabel('range (km)'); xlabel('Doppler vel (m/s)'); end
end
c=colorbar;  set(c,'position',[.925 .1 .01 .2])
%

%%
%print to figure
i=find(filein=='.');

figure(1); hold on
set(gcf,'paperposition',[.25 .15 10 10])
print('-djpeg','-r300',[filein(1:i(1)-1) '_' date '.jpg'])


%%
return

