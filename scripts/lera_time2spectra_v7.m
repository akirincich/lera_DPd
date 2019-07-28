function [SpecHead,data,mtime]=lera_time2spectra_v(incoming_ts_file_dir,incoming_spectra_file_dir,filein,RC,HEAD,CONST)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lera_time2spectra_v.m
%
%   loads a lera timeseries file and processes using the PC structure as
%   well as the information in CONST,HEAD,RC from the processing set up
%   as well as the radar setup file.
% 
%  Variables needed:
%       base_dir
%       filein
%       and the processing constants that will be used:
%           PC.MaxRangeKM=100;        % the maximum ranges, in distance to be kept in spectra, from >0 to Max
%           PC.SpecOverlapPct=50;     %for WOSA, style increased ensembles into CSE
%           PC.SpeclengthPnts=2048;   %defines the number of points in the spectra
%
%
%   This function outputs data to be saved in the matlab file with the filename:
%        CSE_xxx_YYYYyrdhhmm.mat
%
%    where xxx is the site name and yrd is the yearday.
%    This file has
%     mtime    the file start time
%     RC       the radar header information
%     SpecHead the spectrum/data information
%     data     the ansemble averaged auto and cross_spectra from all antennas
%
%   Note that the actual saving of the file is done in the master script
%
%   While the spectra might be estimated using, say 2048, spectral points,
%   only the middle 1024 points are kept moving forward to save space.
%   This is done as the lera chirp is ~1/3 s pushing the bragg region
%   towards the center of the Doppler spectrum, and leaving significant
%   extra space on the outer regions of the spectra that are not required
%   for further analysis
%
%   The output of this function can go right into spectra2radialmetrics steps.
%
%   v2 implemented for realtime processing
%
%   v3 
%   fixes and better explanations. in use for offline, need to
%   transistion.
%
%   v4    8/1/2018
%   moved to more flexible spectral overlapping to allow for overlaps
%   different from 50%, adjusted fft loop to handle this.
%   moved SpecHead.segments to just have the starting index of each ensemble
%   member
%
%  v5    9/2018  changes
%  - changes from the hanning to the Blackman-harris window for improved
%  noise suppression at the expense of a slightly smeared bragg region
%  - also implements hitch finding code to see spikes in the chirp data,
%  and if present, redo the range fft with linear interpolation over hitch
%  location (can be at a different range in each radar, but fairly constant within
%  a file)
%
%    v6  11/2018 changes to streamline lera processing and comment method
%
%    v7  5/2019,  Correcting the velocity issue:  upsweep radars need to flip
%          as the range power will map into the negitive range bin.  Also
%          need to adjust I and Q to match...  This jives with Pierre's
%          recent work on the Taiwanese system...
%           -also fixing a new issue in the hitch detector to allow hitches
%           that are near the close range boundary to be clipped.
%
%        added normalization by rms power at range bin 1 doppler shift ~0
%
% Anthony Kirincich
%   WHOI-PO
%   akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% % %clear
% % %%
% % %%%%% set directories for in script testing %%%%%%%%%%%%%
% % base_dir='~/Matlab/working/LERA/'
% 
% %data_dir
% data_dir=[base_dir 'lera_data/'];
% %script dir
% script_dir=[base_dir 'lera_DP/'];
% %config_dir
% config_dir=[base_dir 'lera_config/'];
% 
% addpath(script_dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%
%unpack processing constants
PC=CONST.PC;
% add to processing constants
% PC.MaxRangeKM=100;        % the maximum ranges, in distance to be kept in spectra, from >0 to Max
% PC.SpecOverlapPct=50;     %for WOSA, style increased ensembles into CSE
% PC.SpeclengthPnts=2048;   %defines the number of points in the spectra

PC.spec_crop=PC.SpeclengthPnts*.25+1:PC.SpeclengthPnts*.75; 

%use to convert from the A2D counts, shifted by the shift bits into a
%voltage
PC.TSunitref=RC.gain/(2^(RC.BIT_LENGTH-RC.SHIFT));    

%PC.TSunitref=10/(2^17);     %use to transform the AD output from counts to volts as 10V = 2^17 counts or 131072
%                         x(in v)/10v = obs_counts/2^17   or x(in v) =
%                         10v*obs_counts/2^17, this is true volts into the receiver
%                         if the receiver is on high gain, max of 10 v
%%% generalize
% switch HEAD.SiteName  %site_name
%     case 'NWTP'
%         %for MK2 receiver
%     PC.TSunitref=RC.gain/(2^(23-RC.SHIFT));     %use to transform the AD output from counts to volts as 10V = 2^17 counts or 131072
%     case 'LPWR'
%         %for MK3 receiver
%     PC.TSunitref=RC.gain/(2^(31-RC.SHIFT));     %use to transform the AD output from counts to volts as 10V = 2^17 counts or 131072
%     case 'HBSR'
%         %for MK3 receiver
%     PC.TSunitref=RC.gain/(2^(31-RC.SHIFT));     %use to transform the AD output from counts to volts as 10V = 2^17 counts or 131072
% end

%%

%%%%% get started, load data file, 
%%% protect RC as the data file might have an older version
RCc=RC;
%%% switch to data dir
eval(['cd ' incoming_ts_file_dir]);
%%% load the file
load(filein);
RC=RCc;  %reload the radar_header version of the RC

%%%% extra if loading a wera ts file
% timechirp=wera; clear wera werac wera1 sdata map data
% mtime=0; file_chirps=2048;

%%%%%%%%%%%%%%%
%flip IQ channels...
timechirp=timechirp([2 1],:,:,:);
swpdir='up';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%build a new header for the spectral output:
SpecHead=[];
SpecHead.SiteName=RC.SiteName;
SpecHead.StartTimeStamp=mtime;          %start time of the dtq file, from the filename
SpecHead.nchirps=file_chirps;           %# of chirps captured in the TS file, could be variable
    
tcsize=size(timechirp);                      %size of the incoming TS data
SpecHead.SamplesPerChirp=tcsize(2);          %#samples per chirp from compressed TS data
SpecHead.RangeMaxKm=RC.c*SpecHead.SamplesPerChirp/4/(RC.BW)./1000;     %max range in km, governed by samples/chirp and BW
SpecHead.RangeResKm=SpecHead.RangeMaxKm*2/SpecHead.SamplesPerChirp;    %the min range resolution, really c/2/BW
%%%%  Where does the '1/4' come from in RangeMax?
%%%%  as done here, this might be a backwards way to establish RangeRes and RangeMax
%%%%  RangeRes is ironclad as c/2/BW, 
%%%%  RangeMax is really a fn of the fft length with is fn of the sample
%%%%  rate AND the rangeres. Thus, the '1/4' as written here comes from the
%%%%  '1/2' for the rangeres * a 1/2 for the sample to fft length.
%%%%      rangemax is better defined as    rangeres*samplesPerChirp/2
SpecHead.FullRangeKm=(-SpecHead.SamplesPerChirp/2:SpecHead.SamplesPerChirp/2-1).*SpecHead.RangeResKm ;  %array of the ranges (neg and positive) for the range FFT

% %if upsweep, reverse range coordinate to account for sweep direction
% %difference
% if strcmp(swpdir,'up')==1;
%     SpecHead.FullRangeKm=-SpecHead.FullRangeKm;
% end

SpecHead.PC=PC;  %copy processing constants  into this header

%work to isolate the range cells that we'll move forward
%just positive range cells:
irange=find((SpecHead.FullRangeKm)<=SpecHead.PC.MaxRangeKM & SpecHead.FullRangeKm>0);
%both pos and negative range cells.
%irange=find(abs(SpecHead.FullRangeKm)<=SpecHead.PC.MaxRangeKM);

SpecHead.RangeKm=SpecHead.FullRangeKm(irange);
SpecHead.rangecell=1:length(irange);

%establish frequency, period
SpecHead.Fc=RC.Fc;                %center frequency of radar transmission
if SpecHead.Fc < 1e6;  %this is not in Hz, assume its in MHz and mult by 1e6
    SpecHead.Fc=SpecHead.Fc.*1e6;   % now in Hz
end
SpecHead.Tr=RC.chirp;              %chirp period s
SpecHead.fr=1./RC.chirp;            %the chirp rate in Hz
SpecHead.fDmax=.5*SpecHead.fr;    %max dop frequency resolution, for the spectral estimate itself

% more frequencies
%  DopplerResolution = 1.0 / AveragingPeriod
SpecHead.delta_f=2*SpecHead.fDmax/(SpecHead.PC.SpeclengthPnts);     % correct ..... 
SpecHead.doppler_freq=(-SpecHead.fDmax+SpecHead.delta_f:SpecHead.delta_f:SpecHead.fDmax);  %correct......

%%% find what the bragg frequency is
SpecHead.L=(RC.c/SpecHead.Fc)/2;    %lamda/2
SpecHead.v_p=sqrt(9.81*SpecHead.L/2/pi);
SpecHead.FBragg=2*SpecHead.v_p*SpecHead.Fc/RC.c;
[s,i]=sort(abs(abs(SpecHead.doppler_freq)-SpecHead.FBragg));
SpecHead.iFBragg=i(1:2);

%get doppler_vel
SpecHead.doppler_vel=SpecHead.doppler_freq*RC.c/SpecHead.Fc/2;
%c_vel=doppler_vel-[-v_p*ones(1,(spectra_length-2)/2) 0 v_p*ones(1,(spectra_length-2)/2)];
SpecHead.c_vel=SpecHead.doppler_vel-[-SpecHead.v_p*ones(1,(SpecHead.PC.SpeclengthPnts)/2) SpecHead.v_p*ones(1,(SpecHead.PC.SpeclengthPnts)/2)];

%%% narrative for what all this does %%%%%
%VelocityResolution = c * DopplerResolution / TransmitCenterFreqMHz./1e6
%VelocityMaximum = c * DopplerMaximum / TransmitCenterFreqMHz./1e6
%DopplerMaximum = DopplerResolution * NumberDopplerCells / 2

%establish normpoint, the zero range, zero freq point to normalize by
normpt=[find(SpecHead.RangeKm==min(abs(SpecHead.RangeKm))) find(SpecHead.doppler_freq==0)];   


%%
%%%%%%%%%% move forward and process the TS data to the range and spectral estimates
%make the range data array
range_data=zeros(SpecHead.SamplesPerChirp,SpecHead.nchirps,RC.NANT);
%%% if investigating the position of the chirp max
%outto_chirpmax=[];

% %lera often has a hitch in the chirps at a fixed position, when this
% %occurs interp through a small gap reduces noise by 10db, When it is not
% %present the interp seems to have little effect. 
% v5 of this script implements a redo of the range fft if a hitch is found
% and replots the difference.
%
%  This is a patch for fixing something about the radar itself.
dysave=nan.*ones(SpecHead.SamplesPerChirp,file_chirps);  %create variable for saving the mean abs(chirp) to find a hitch
hitch_factor=4;
hitch_factor_span=3;

%%% loop over chirps and compute range spectra
for ichirp=1:file_chirps
    %isolate the i and q parts for this chirp, place into units of volts.
    wc=double(timechirp(:,:,:,ichirp)).*SpecHead.PC.TSunitref;   
    %form the complex variable
    y=double(squeeze(wc(1,:,:))) + sqrt(-1).*double(squeeze(wc(2,:,:))) ;  %combine into a complex #
    dysave(2:end,ichirp)=mean([abs(diff(real(y))) abs(diff(imag(y)))],2);     % get data for hitch finder
    %do fft
%    Y=fft(y.*(hanning(SpecHead.SamplesPerChirp)*ones(1,RC.NANT)));   %do the fft with a hanning window .  %better side lobe suppression.
    Y=fft(y.*(blackmanharris(SpecHead.SamplesPerChirp)*ones(1,RC.NANT)));   %do the fft with a hanning window .  %better side lobe suppression.
%    %shift to put the 0 frequency in the middle of the array 
    Y=fftshift(Y,1);  %in agreement with SpecHead.FullRangeKm
    
    %if upsweep, reverse range coordinate to account for sweep direction
    %difference
    if strcmp(swpdir,'up')==1;
        Y=flipud(Y);
    end
    
    %%% place this result into new array of 
    range_data(:,ichirp,:)=Y;

%     %correct the variance? out side of the big jump
%     fs=1./(SpecHead.Tr/SpecHead.SamplesPerChirp)  %freq increment 
%     parPy=(var(real(y(20:end-20,:))).*((SpecHead.SamplesPerChirp-1)./SpecHead.SamplesPerChirp))./abs(sum(real(Y).*fs))  
    
%%%% make a plot of the range result for this chirp 
if CONST.goplot(2)==1;
    if ichirp==1 | ichirp==file_chirps;  %only plot the first and last
    figure(10); clf;
    subplot(211);  plot(real(y)); hg; plot(imag(y));
    a=ceil(SpecHead.SamplesPerChirp.*[.2 .5]);
    axis([1 SpecHead.SamplesPerChirp  nanmax(nanmax(abs(real(y(a(1):a(2),:)))))*[-3 3]])
    title('Timeseries for this chirp');
    subplot(223);
    plot(SpecHead.FullRangeKm,20.*real(log10(real(Y))),'k'); hg; 
    axis([[-1 1].*max(SpecHead.FullRangeKm) -90 20])
    title('I channel range dependent power'); xlabel('km')
    subplot(224)
    plot(SpecHead.FullRangeKm,20.*real(log10(imag(Y))),'r'); hg;
    title('Q channel range dependent power'); xlabel('km')
    axis([[-1 1].*max(SpecHead.FullRangeKm) -90 20])
    pause(.05)
    end
end
% %check on the chirp timing
% rout=abs(real(y))./( ones(length(y),1)*max(abs(real(y))) );
% iout=abs(imag(y))./( ones(length(y),1)*max(abs(imag(y))) );
% 
% i=find(mean(rout,2)==max(mean(rout,2)));
% j=find(mean(iout,2)==max(mean(iout,2)));
% outto_chirpmax(ichirp,:)=[i(1) j(1) max(mean(rout,2)) max(mean(iout,2))];

end  %for ichirp

%%

if CONST.gohitch==1
%%% examine data for hitches:
y=nanmean(dysave,2);
ia=floor(SpecHead.SamplesPerChirp*.15):SpecHead.SamplesPerChirp-floor(SpecHead.SamplesPerChirp*.10);
me=nanmean(y(ia)); ma=nanmax(y(ia));

if CONST.goplot(2)==1; figure(20); clf; plot(y); hg; title('hitch finder'); end

if ma>me.*hitch_factor; % a persistent hitch in the chirp is present.  redo this work with
    
i=find(y(ia)==ma);
%    gap=ia([i-hitch_factor_span i+hitch_factor_span]);
    gap=[ia(i)-hitch_factor_span ia(i)+hitch_factor_span];  %fix on 5/2019
    igap=[(1:gap(1)-1) (gap(2)+1:SpecHead.SamplesPerChirp)];
    comment=sprintf('A hitch between chirp indices %2.0f and %2.0f was interped over',gap); 
    disp(comment)

%make the range data array..again with the hitch interped
range_data=zeros(SpecHead.SamplesPerChirp,SpecHead.nchirps,RC.NANT);

%%% loop over chirps and compute range spectra
for ichirp=1:file_chirps
    %isolate the i and q parts for this chirp, place into units of volts.
    wc=double(timechirp(:,:,:,ichirp)).*SpecHead.PC.TSunitref;   
    %form the complex variable
    y=double(squeeze(wc(1,:,:))) + sqrt(-1).*double(squeeze(wc(2,:,:))) ;  %combine into a complex #
    %interp over gap area of known bad data...  
    y=interp1(igap',real(y(igap,:)),1:SpecHead.SamplesPerChirp) + sqrt(-1).*interp1(igap',imag(y(igap,:)),1:SpecHead.SamplesPerChirp) ;
    %do fft
%    Y=fft(y.*(hanning(SpecHead.SamplesPerChirp)*ones(1,RC.NANT)));   %do the fft with a hanning window .  %better side lobe suppression.
    Y=fft(y.*(blackmanharris(SpecHead.SamplesPerChirp)*ones(1,RC.NANT)));   %do the fft with a hanning window .  %better side lobe suppression.    %    %shift to put the 0 frequency in the middle of the array  %the wrong way to shift?
   Y=fftshift(Y,1);  %in agreement with SpecHead.FullRangeKm
     %if upsweep, reverse range coordinate to account for sweep direction
    %difference
    if strcmp(swpdir,'up')==1;
        Y=flipud(Y);
    end
    
    %%% place this result into new array of 
    range_data(:,ichirp,:)=Y;

%%%% make a plot of the range result for this chirp 
if CONST.goplot(2)==1;
    if ichirp==1 | ichirp==file_chirps;  %only plot the first and last
    figure(10); clf;
    subplot(211);  plot(real(y)); hg; plot(imag(y));
    a=ceil(SpecHead.SamplesPerChirp.*[.2 .5]);
    axis([1 SpecHead.SamplesPerChirp  nanmax(nanmax(abs(real(y(a(1):a(2),:)))))*[-3 3]])
    title('Timeseries for this chirp');
    subplot(223);
    plot(SpecHead.FullRangeKm,20.*real(log10(real(Y))),'k'); hg; 
    axis([[-1 1].*max(SpecHead.FullRangeKm) -90 20])
    title('I channel range dependent power'); xlabel('km')
    subplot(224)
    plot(SpecHead.FullRangeKm,20.*real(log10(imag(Y))),'r'); hg;
    title('Q channel range dependent power'); xlabel('km')
    axis([[-1 1].*max(SpecHead.FullRangeKm) -90 20])
    pause(.05)
    end
end

end %for ichirp
end %if hitch is present
end  %if gohitch

%%% test case for proving fft and fftshift, following matlab docs
%  
% fs = 100;               % sampling frequency
% t = 0:(1/fs):(10-1/fs); % time vector
% S1 = cos(2*pi*15*t);
% S2 = cos(2*pi*30*t);
% n = length(S1);
% A = [S1; S2];
% X = fft(A,[],2);
% f = (0:n-1)*(fs/n);     % frequency range
% power = abs(X).^2/n;    % power
% plot(f,power(1,:),f,power(2,:))
% 
% Y = fftshift(X,2);
% fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
% powershift = abs(Y).^2/n;     % zero-centered power
% plot(fshift,powershift(1,:),fshift,powershift(2,:))

%%%%%%%%

%%% make a plot of the power of the I and Q in all chirps for ant 1
if CONST.goplot(2)==1;
    figure(11); clf;
    subplot(211);  pcolor(20.*real(log10(real(range_data(:,:,5))))); shading flat; colorbar; hg;
    subplot(212);  pcolor(20.*real(log10(imag(range_data(:,:,5))))); shading flat; colorbar;  hg;

%     figure(11); clf;
%     plot(outto_chirpmax(:,1)); hg;
%     plot(outto_chirpmax(:,2))
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% prep for second FFT to spectra
%cut range data down to +/- range of interest.
range_data_cut=range_data(irange,:,:);
%figure out how many spectra will be averaged together to make the ensemble estimate
%%% for 50% overlap
%SpecHead.segments=1:SpecHead.PC.SpeclengthPnts*((100-SpecHead.PC.SpecOverlapPct)/100):SpecHead.nchirps;
%SpecHead.nsegs=length(SpecHead.segments)-2;

%%% for random overlap
SpecHead.segments=floor(1:SpecHead.PC.SpeclengthPnts*((100-SpecHead.PC.SpecOverlapPct)/100):SpecHead.nchirps-SpecHead.PC.SpeclengthPnts+1);
SpecHead.nsegs=length(SpecHead.segments);


%%% only proceed if there are enough samples to make one ensemble member
if SpecHead.nsegs > 1

%make outgoing variable
saveY=nan.*ones(length(irange),SpecHead.PC.SpeclengthPnts,RC.NANT,SpecHead.nsegs);
aa=nan*ones(1,8); 

%%% do second FFT for each ant and segment, separately,
for ii=1:RC.NANT
    for iii=1:SpecHead.nsegs;
        b=SpecHead.segments(iii):SpecHead.segments(iii)+SpecHead.PC.SpeclengthPnts-1;   %adjusted v4 for flexibility
        y=range_data_cut(:,b,ii);
        Y=fft(y.*(ones(length(irange),1)*(hanning(SpecHead.PC.SpeclengthPnts)')),[],2);   %do the fft with a hamming window   
        Y=fftshift(Y,2);  %shift to put the 0 frequency in the middle.
        %       saveY(:,:,ii,iii)=Y;
   %try normalizing each antenna by the rms power in a range around the center point  
        nabs=sqrt(nanmean(abs(Y(normpt(1),[normpt(2)-floor(SpecHead.PC.SpeclengthPnts/500):normpt(2)+floor(SpecHead.PC.SpeclengthPnts/500)])).^2));
        mag=abs(Y(:,:));
        pha=atan2(imag(Y),real(Y));
        saveY(:,:,ii,iii)=(mag./nabs).*exp(sqrt(-1).*pha);          
    end
    
  %%Make plots of the spectra as they are created
 if CONST.goplot(2)==1
    figure(2); clf;
    for i=1:SpecHead.nsegs;
        subplot(SpecHead.nsegs+1,1,i);
        pcolor(20*log10(abs(squeeze(saveY(:,:,ii,i))))); shading flat; colorbar; hg; caxis([-100 20]); hg
        plot([1; 1]*SpecHead.iFBragg,SpecHead.rangecell([1 end])'*[1 1],'k')
        title(['Antenna ' num2str(ii) '; Segment ' num2str(i)])
    end
        subplot(SpecHead.nsegs+1,1,SpecHead.nsegs+1);
        pcolor(20*log10(squeeze(nanmean(abs(saveY(:,:,ii,:)),4)))); shading flat; colorbar; hg; caxis([-100 20]); hg;
        plot([1; 1]*SpecHead.iFBragg,SpecHead.rangecell([1 end])'*[1 1],'k')
        title(['Antenna ' num2str(ii) '; average '])

        i=find((SpecHead.c_vel-floor(SpecHead.c_vel))<0.02);
        set(gca,'xtick',i,'xticklabel',num2str(round(SpecHead.c_vel(i))'))
        
        aa(ii)=nanmean(nanmean(20*log10(squeeze(nanmean(abs(saveY(:,:,ii,:)),4)))));
  pause(.5)
 end  %if goplot
end

%%

%%%%%% create the ensemble-averaged auto- and cross-spectral estimates
%cut the outer parts of the spectra for space using PC.spec_crop
SpecHead.c_velc=SpecHead.c_vel(SpecHead.PC.spec_crop);
SpecHead.doppler_velc=SpecHead.doppler_vel(SpecHead.PC.spec_crop);
SpecHead.doppler_freqc=SpecHead.doppler_freq(SpecHead.PC.spec_crop);
SpecHead.iFBraggc=SpecHead.iFBragg-SpecHead.PC.spec_crop(1)+1;

%resave Y as the smaller size
Y=saveY(:,SpecHead.PC.spec_crop,:,:);
%%%%%% drop this into a form similar COS with the 
data=[];

for ii=1:RC.NANT;
    %make the autospectra for this antenna
% %    eval(['data.a' num2str(ii) num2str(ii) '=sum(abs(saveY(:,:,ii,:)).^2,4)./SpecHead.nsegs;']) 
%      eval(['data.a' num2str(ii) num2str(ii) '=sum(abs(Y(:,:,ii,:)).^2,4)./SpecHead.nsegs;'] );
      eval(['data.a' num2str(ii) num2str(ii) '=sum( [Y(:,:,ii,:).*conj(Y(:,:,ii,:))],4)./SpecHead.nsegs;'] );
    %make the cross-spectra for this antenna against all the higher antennas
    for jj=ii+1:RC.NANT;
%        eval(['data.a' num2str(ii) num2str(jj) '=sum(saveY(:,:,ii,:).*conj(saveY(:,:,jj,:)),4)./SpecHead.nsegs;']) 
        eval(['data.a' num2str(ii) num2str(jj) '=sum(Y(:,:,ii,:).*conj(Y(:,:,jj,:)),4)./SpecHead.nsegs;'] );
    end
end


else;  data=[];  
    disp('this file is short on chirps, returning and empty array for data')
end  %SpecHead.nsegs > 1




return








