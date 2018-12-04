function [R, HEAD]=HFR_DP_lera_music(data,FOregi,HEAD,SpecHead,patt,CONST)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [R, HEAD]=HFR_DP_lera_music(Data,FOregi,HEAD,SpecHead,patt,CONST)
%
%  follows MUSIC algorithm (Schmidt,1986) to compute and
%    evaluate the DOA function for each spectral point in FOregi for up to
%    Nmax solutions.  For each MUSIC result, the DOA
%      function is evalulated to find the significant peaks using
%      pksfinder and a super-complex set of conditional statments to
%      account for a variety of measured pattern types and azimuthal extents
%      (see below)
%
% INPUT:  Should be self-explanatory given the prepatory work by the
% spectra2radialmetric script that calls this function.
%
%outputs:
%
%   R  --  A matrix of radial metric output.  This output
%          follows COS documentation for column/field order for lack
%          of a better option. Outgoing data starts with the range
%          (column 6) of COS radial metrics file. The rest (lon lat
%          u v flag) are added later.
%
%         The number rows of R is variable and depends on the number of
%         FOregi and the percent of duel angle solutions returned.
%
%%%% put outgoing data into an array similar to radial metrics, starting
%%%%  add the rest (lon lat u v flag) after all FOLs are processed.
%%%%  (lon lat u v flag) will be added to the front end later
% 1 range
% 2 bearing
% 3 vel
% 4 direction
% 5  rangecell
% 6  dopcell
% 7  angselect (which solution if isave out is given )
% 8  which peak of this solution is given
% 9 musicpow(v)
% 10 music_pkwidth  %these are really wide for low angselect numbers, (in seasonde, this would mean a bad pattern?)
% 11 musicDOASpeak
% 12 spectraA3(snr)
% 13-15 musicEigenRatio musicpowerRatio musicoffRatio
% 16-22 which solutions were viable (1 yes 0 no, for 1-7 (need better description/name)
% 23-30 Eigenval (1-8)
% 31   DF_flag (see above for explaination
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
% Versions:
% . created 8/2017
%
%   edited 4/2018
%    fixed max of single angle solution to be single peak,
%   i.e. peak finding is not done on the single angle case and this case is
%   always valid.
%
%
%   edited 5/2018 to allow for two types of emitter number methods to be
%   used.  These are based on switches in CONST and either can be used
%
%   CONST.which_Ns_meth='music_param';
%   CONST.which_Ns_meth='music_highest';
%
%   also this requires the constant   CONST.doa_peak_thresh   to be set.
%
%   This defines the peak threshold to use for the music doa peak finder.  Note
%   that this is done in 10*log10(MS) space, making .05-.5 reasonable thresholds
%   that scale with the values to allow small peaks to exist.
%
%  v2  11/2018 clean up and consolidation for similarity with
%  HFR_DP_lera_DF_v?.m
%
%    Anthony Kirincich
%    WHOI-PO
%    akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%% for debuging
%MSS_save=[];


%start with arrays for all possible angle solutions R1-RNmax
% condense and combine after finished
for ii=1:CONST.Nmax
    eval(['R' num2str(ii) '=nan.*ones(length(FOregi),31);'])   %output array for sing or Duel1
    %R2=nan.*ones(length(FOregi),29);   %output for Duel2
end

%%% prep for calculations
%get the mean powerspectra
gains=[];
for i=1:CONST.N
    eval(['gains(:,:,i)=(data.a' num2str(i) num2str(i) '); '])
end
%%%%% establish the noise floor
%%% find the one with the lowest noise floor and use this as the noise floor
noise_floor=nanmean(squeeze(gains(end,:,:)));
i=find(noise_floor==min(noise_floor));
noise_floor=noise_floor(i(1));
gain3=gains(:,:,i(1)); clear gains;

%%%% note that HFR_DP_lera_DF_v? uses the slightly different method of
%%%% looking at the mean noise floor across all ants.
%%%%  This is used as
% %%% method that uses log10 dB scale to establish the noise floor.
% gain3=10*log10(nanmean(gains,3)) + (-40 - 5.8); clear gains;
% % use the mean power in the outer edge to get the noise floor
% noise_floor=nanmean(gain3(end,[1:CONST.imageFOL_user_param end-CONST.imageFOL_user_param:end]));


%extract constants.
mP=HEAD.Musicparams123(1:3);
peak_thresh=CONST.doa_peak_thresh;
dangles=nanmedian(abs(diff(patt.angles)));

%%
%%%   loop over all rows of FOregi and compute DOA and MUSIC signal qualities for each.
%do direction finding for each good FOregi
for i4=1:1:length(FOregi)
    %%
    %isolate point of interest for this loop
    bp=FOregi(i4,:);
    
    %%%% start music script %%%%
    %sends the covariance, the APM, and theta, and the max number of signals
    %returns the peaks of each DOA, up to N, after evaluation for peaks on the boundaries.
    %(follow existing music code)
    
    %%% (1) make covariance function,
    %loop over all antennas to fill out A
    C=nan.*ones(CONST.N,CONST.N);
    
    for ii=1:CONST.N
        eval(['C(ii,ii)=data.a' num2str(ii) num2str(ii)   '(bp(1),bp(2));'])
        for jj=ii+1:CONST.N
            eval(['C(ii,jj)=data.a' num2str(ii) num2str(jj)   '(bp(1),bp(2));'])
            eval(['C(jj,ii)=conj(data.a' num2str(ii) num2str(jj)   '(bp(1),bp(2)));'])
        end
    end
    
    if isempty(find(isnan(C(:))==1))==1  %fix to ensure that data in C is good
        
        %%%% (2) Estimate DOA function
        %%%% equation to solve is c = a*s*a*T
        [Q ,D]=eig(C); %Compute eigendecomposition of covariance matrix
        [D,I]=sort(diag(D),1,'descend'); %Find r largest eigenvalues
        Q=Q (:,I); %Sort the eigenvectors to put signal eigenvectors first
        
        MS=[];
        %compute music for N=1:max(N)-N/4 DOAs,  upper solutions get very noisy
        for j=1:CONST.Nmax;
            %Qs=Q (:,1:r); %Get the signal eigenvectors, not actually needed
            Qn=Q(:,j+1:CONST.N); %Get the noise eigenvectors
            music_spectrum=[];
            for i=1:length(patt.A)
                %%%find the value for each direction.
                %            music_spectrum(i)=diag( [ (patt.A(:,i)'*patt.A(:,i))/(patt.A(:,i)'*Qn*Qn'*patt.A(:,i)) ]);
                %         %%% Compute the DOA function (Schmidt 1986, equation 6)
                music_spectrum(i)=( [ (patt.A(:,i)'*patt.A(:,i))/(patt.A(:,i)'*Qn*Qn'*patt.A(:,i)) ]);
                %        DOA(j,i) = 1/(A(:,i)'*(Qn*Qn')*A(:,i)); %does same above, just different by a factor of ~9 for an ideal ant
            end
            MS(j,:)=real(music_spectrum);
        end
        
        %%% (3) get the peaks of all MS, As I'm not often using a, ignore all
        %%% peaks on the boundary, but be sure to place the boundaries on the
        %%% back,land side of the pattern as arranged around the antenna array
        ipks=cell(CONST.Nmax,1);
        igood=nan*ones(1,CONST.Nmax);
        %assume that max of 1 emitter solution is the only peak for this solution
        ipks(1)={find(MS(1,:)==max(MS(1,:)))}; igood(1)=1; % is the only peak,
        %thus single angle solution is always valid... move on to the others
        for j=2:CONST.Nmax;
            [pks,dzdt]=pksfinder(10*log10(MS(j,:)),peak_thresh);
            %focus on the upper, isolated peaks
            ipks(j)={find(diff(dzdt)==-2)+1};
            %%%  which MS have the correct number of peaks?
            if length(ipks{j})==j; igood(j)=1; end
        end
        
        if CONST.goplot(2)==1
            %%% plot this result so far %%%
            figure(4); clf;  colors='krbgc';
            for j=1:CONST.Nmax;
                %            plot(patt.angles,(MS(j,:)'),[colors(j) '.'],'markersize',4); hg;  %plot in regular space
                %            plot(patt.angles(ipks{j}),(MS(j,ipks{j})),'ko');
                plot(patt.angles,10*log10(MS(j,:)'),[colors(j) '.'],'markersize',4); hg;  %plot in log space
                plot(patt.angles(ipks{j}),10*log10(MS(j,ipks{j})),'ko');
                
                if igood(j)==1;
                    %                plot(patt.angles,(MS(j,:)'),[colors(j),'.'],'markersize',12);  %plot in regular space
                    plot(patt.angles,10*log10(MS(j,:)'),[colors(j),'.'],'markersize',12); %plot in log space
                end
            end
            plot(patt.angles,5*ones(size(patt.angles)),'k')
            axis([0 360 0 25])
            title(num2str(bp))
            pause; %(.2)
        end
        %%%%%%%%%%% end music script for this FOregi %%%%%%%%
        
        %%
        %%%%%%%%%%%  start evaluation script. %%%%%%%%
        
        %%% (4) evaluate the DOA function to assess how many scatterers there are
        %%% steps to id which MS is the correct one
        
        %%% (A) cut DOAs to MS where # peaks matches predicted, (using igood)
        %%%     as we shouldn't trust a 3 scatter solution that only returns 2
        %%%     scatterers
        isave=find(igood==1);
        pks=ipks(isave);
        mMS=MS(isave,:);
        
        %%%%%% continue to calculate power, and output variables %%%%
        %%%% in schmidt notation %%%
        %  C is the covariance matrix
        %   Co is the noise covariance matrix but unknown...assume it is the identity matrix
        %    lambda_min is the noise eigenvalue
        %Co=[1 0 0; 0 1 0; 0 0 1]; %assume the noise covariance matrix is the identity matrix
        Co=zeros(CONST.N); for i=1:CONST.N; Co(i,i)=1; end
        
        if isempty(isave)==0;   %only proceed if there is one, viable solution
            %%% loop over all isave and export the Schmidt estimated power.
            clear P1 P2 P3 P4 P5 P6
            for ii=1:length(isave)
                %use found peak directions to calculate A and P (in schmidt) or S (in lipa et al)
                itheta=pks{ii}';
                %     %sort itheta by the DOA magnitude
                [s,i]=sort(mMS(ii,itheta),'descend');
                itheta=itheta(i);
                A1=patt.A(:,itheta);
                eval(['P' num2str(ii) '=((A1''*A1)^-1)*A1''*(C- D(isave(end)+1)*Co)*A1*((A1''*A1)^-1);'])  %back out the signal power
                %put the reordered itheta back into the pks data base...so that order of
                %pks at the end matches the doa mag.
                pks{ii}=itheta;
                
            end
            
            
            %%%%%%%%%%%
            %%% form the exporting structure, similar for all parts below
            M=[];
            M.Eigen=[abs(D)];   % we should give all the eigenvalues every time
            %find mean snr of all antennas
            M.snrant=[mean(diag(C)./noise_floor) ];
            M.isaveout=zeros(1,CONST.N-1);
            M.isaveout(isave)=1;
            M.param=[nan nan nan];   %starter set for these, will change ony if 2+ solutions
            M.crit=[nan nan nan];
            M.DF_flag=nan;
            jj=1; %start with P1 as the winner and test
            
            %%%%%%%%%%%%%%%%
            if length(isave)>1;  %2+ matches, compare 1 to 2 and winner to 3rd, serially
                
                %%% revise flag for multiple angle solutions
                M.DF_flag=0;
                
                %%% Start the evaluation part of this work if there are multiple good
                %%% solutions.
                if strcmp(CONST.which_Ns_meth,'music_highest')==1;  % pick the highest solution that has the right num peaks
                    jj=length(isave);   %done, let's get outa here
                    
                elseif strcmp(CONST.which_Ns_meth,'music_param')==1;  %follow lipa to use params to eval solutions.
                    
                    %compare the first to the second and then the winner to the third, etc.
                    for ii=1:length(isave)-1
                        %prep
                        eval(['p=P' num2str(ii+1) ';']);  % make p the signal matrix for the higher solution
                        d=sort(abs(diag(p)),'descend');   %pull out and organize the diagonals
                        for i=1:length(p);  p(i,i)=nan + sqrt(-1).*nan;  end  %strip the diags away for the calculation
                        %get music parameters
                        M.param=[abs(D(isave(jj)))./abs(D(isave(ii+1))) ...   %the ratio of the eigenvalues
                            mean( d(1:end-1)./d(2:end)  )      ...       %the ratio of the signal powers from P*
                            abs(prod(d)) / abs(prod(p(:),'omitnan'))];   %the ratio of the diags to off-diags
                        M.crit=[0 0 0];
                        if M.param(1) < mP(1);        M.crit(1)=1;    end
                        if M.param(2) < mP(2);        M.crit(2)=1;    end
                        if M.param(3) > mP(3);        M.crit(3)=1;    end
                        
                        %%%%   make the flag of the DF parameters
                        %%% The old flag  was M.DF_flag=(7-(M.crit*[1 2 4]')); for each iteration of the comparison.
                        
                        %%% since the successful number of peak info is elsewhere, DF_flag only
                        %%% needs to show why a single/duel (for example) was made
                        %combine all flags together (each get 2 digits)
                        % if the flag is nan, only one solution existed
                        % if the flag is >0,  reject the duel angle solution in favor of the
                        % single or first
                        M.DF_flag=M.DF_flag+(7-(M.crit*[1 2 4]')).*10^(ii-1);
                        if (7-(M.crit*[1 2 4]'))==0; jj=ii+1;   %use the duel or higher angle solution,
                            % follow the DF_flag way of doing this
                        end
                    end  %loop over pairs of solutions in order.
                end  % if strcmp
            end  %if length (isave) greater than 1.
            
            %get the final itheta
            itheta=pks{jj}';
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% loop over the number of peaks for this solution and export the
            %%% results into one line of RM file for each solution
            for iii=1:length(itheta)
                
                %%%%%%%% gather peak specific data  %%%%%%%%%
                M.angles=[patt.angles(itheta(iii))];
                eval(['M.power=P' num2str(jj) '(iii,iii);']);
                M.DOApeak=abs([mMS(jj,itheta(iii)) ]);
                
                %%% find the 3dB width of this peak, this can be tricky if the peak is
                %%% small or there are other peaks around, that are larger than this
                %%% one.
                a=3; %3 dB by default
                if range(real(mMS(jj,:)))<a   %sometimes the whole thing is quite small, adjust the 3db downwards
                    a=max([1 range(real(mMS(jj,:)))./2]);
                end
                j=find(real(mMS(jj,:))>real(mMS(jj,itheta(iii)))-a);
                M.halfpowwidth(1)=abs(length(j)./dangles);
                %%% if jj>1 assume that the other peaks might be interfering with this
                %%% half power width. thus constrain the width to the length at a down
                %%% which is j OR the distance between the peaks if smaller, as a proxy
                %%% for the width.
                if jj>1;  % if more than one angle exists.
                    angle_gaps=itheta-itheta(iii); angle_gaps(find(angle_gaps==0))=nan  ;
                    i=find(abs(angle_gaps)==nanmin(abs(angle_gaps)));
                    if length(j)>abs(angle_gaps(i))
                        M.halfpowwidth(1)=abs(angle_gaps(i(1))./dangles);
                    end
                end
                
                %    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%% put outgoing data into an array similar to the radial metrics format, starting
                %%% add the rest (lon lat u v flag) after all FOLs are processed.
                % 1 range
                % 2 bearing
                % 3 vel
                % 4 direction
                % 5  rangecell
                % 6  dopcell
                % 7  angselect (which solution of isave out is given )
                % 8  which peak of this solution is given
                % 9 musicpow(v)
                % 10 music_pkwidth
                % 11 musicDOASpeak
                % 12 spectraA3(snr)
                % 13-15 musicEigenRatio musicpowerRatio musicoffRatio
                % 16-22 which solutions were viable (1 yes 0 no, for 1-7 (need better description/name)
                % 23-30 Eigenval (1-8)
                % 31   DF_flag (see above for explanation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%% export data section %%%%%%%%%
                % d is  [ range bearing vel direction rangecell dopcell angleselect]
                d=[SpecHead.RangeKm(bp(1))  M.angles SpecHead.c_velc(bp(2))  M.angles+180];    if d(4)>360; d(4)=d(4)-360; end
                r=[d bp jj iii abs(M.power) M.halfpowwidth M.DOApeak M.snrant M.param M.isaveout M.Eigen' M.DF_flag];
                eval(['R' num2str(iii) '(i4,:)=r;'])
                
                %  jj and iii tell you what solution was chosen and what angle of this
                %  solution is given in the line of data.
                %%%%%%%% end export data section %%%%%%%%%
                
                if CONST.goplot(2)==1 & sum(M.isaveout)>1;   M
                    pause
                end
                
            end % for iii length(isave) loop over all peaks in the single solution
        end %if isempty isave
    end  %if isempty C
end  %for i4 loop over all FOregi


%combine all Rs, while getting rid of empty lines (marked by nans)
for ii=1:CONST.Nmax
    eval(['i' num2str(ii) '=find(isnan(R' num2str(ii) '(:,1))==0);'])
end
if CONST.Nmax==3
    R=[R1(i1,:); R2(i2,:); R3(i3,:)];
elseif CONST.Nmax==4
    R=[R1(i1,:); R2(i2,:); R3(i3,:); R4(i4,:)];
elseif CONST.Nmax==5
    R=[R1(i1,:); R2(i2,:); R3(i3,:); R4(i4,:); R5(i5,:)];
end
clear R1 R2 R3 R4 R5

%

%%%%% add this script to the processing steps
HEAD.ProcessingSteps{end+1}=mfilename;

return

