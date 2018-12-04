function [R, HEAD]=HFR_DP_lera_df_vd(data,FOregi,HEAD,SpecHead,patt,CONST)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [R, HEAD]=HFR_DP_lera_df_vd(data,FOregi,HEAD,SpecHead,patt,CONST)
%
%   called from HFR_DP_lera_spectra2radialmetric_process_v
%     this is a functional script to load Hermitian covariance matrix for each FOregi point
%     to consider for direction finding, estimate the number of emitters present
%     and use a DF method to determine directions of arrival.
% 
%    HFR_DP_lera_df_v  has the ability to use music, mle and wsf methods to 
%     determine the doas.
%
%    each method returns a doa function (or the trace), the indices of the
%    estimated emitter locations, and the iteration number (for mle and
%    wsf)
% 
%   DF method results are interrogated to fill out the power, signal
%   peak/width information for each emitter. Each emitter is reported in a
%   line of the returning matrix R along with information about the 
%
%   The estimated number of emitters is found using the statisitical methods 
%    described by Johnson and Dudgen, sec 7.3.5
%
%
% INPUT:  Should be self-explanatory given the prepatory work by the
% spectra2radialmetric script that calls this function.
%
% OUTPUT:
%
%   R  --  A matrix of radial metric output starting with the range
%           of radial metrics file. The rest (lon lat
%          u v flag) are added later.
%
%         The number rows of R is variable and depends on the number of
%         FOregi and number of multiple angle solutions
%
%       R columns given here are:
%           1 range
%           2 bearing
%           3 vel
%           4 direction
%           5  rangecell
%           6  dopcell
%           7  total number of emitters for this FOreg
%           8  which peak of this solution is given
%           9 musicpow(v)
%           10 music_pkwidth  %these are really wide for low angselect numbers, (in seasonde, this would mean a bad pattern?)
%           11 musicDOASpeak
%           12 spectraA3(snr)
%           13-15 musicEigenRatio musicpowerRatio musicoffRatio
%           16-22 which solutions were viable (1 yes 0 no, for 1-7 (need better description/name)
%           23-30 Eigenval (1-8)
%           31   DF_flag (see above for explaination
%
%
%    As a reminder, the N=8 lera rectangular array should be thought of in
%    the following coordinate system... 
%
%        ^y          in math coord from the pos x axis,
%        |      pos bearings
%        o   o
%  _ o   o   o  _ > x
%    o   o   o
%        |        neg bearings
%
%  at nwpt this is a bearing of 147+90T
%
%
% Versions:
%   created 8/2017
%   altered to allow multiple methods, 12/2017
%
%    Anthony Kirincich
%    WHOI-PO
%    akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%start with arrays for all possible angle solutions R1-RNmax
% condense and combine after finished
for ii=1:CONST.Nmax
    eval(['R' num2str(ii) '=nan.*ones(length(FOregi),31);'])   %output array for sing or Duel1
end

%%% prep for calculations %%%%%
%%% establish the noise floor, by looking at the spectral power levels
gains=[];
for i=1:CONST.N
    eval(['gains(:,:,i)=abs(data.a' num2str(i) num2str(i) '); '])
end

%%%%% establish the noise floor %%%%%%%%
%%% find the one with the lowest noise floor and use this as the noise floor
noise_floor=nanmean(squeeze(gains(end,:,:)));
i=find(noise_floor==min(noise_floor));
noise_floor=noise_floor(i(1));
gain3=gains(:,:,i(1)); clear gains;

%%% use an alternative that focus on the mean noise floor of all antennas
%%%,  AND keeps the gain3 and noise_floor in dB, not volts
%%% if using, we need to adjust the M.snrant code below around line 244

% % method that uses log10 dB scale to establish the noise floor.
% gain3=10*log10(nanmean(gains,3)) + (-40 - 5.8); clear gains;
% % use the mean power in the outer edge to get the noise floor
% noise_floor=nanmean(gain3(end,[1:CONST.imageFOL_user_param end-CONST.imageFOL_user_param:end]));


%%% set coefficients for direction_finding, etc.
peak_thresh=CONST.doa_peak_thresh;   %do this in 10*log10(MS) space, making .5 a reasonable threshold?
dangles=nanmedian(abs(diff(patt.angles)));
mP=[];

%%%   loop over all rows of FOregi and compute DOA and MUSIC signal qualities for each.
%%

%%% prep for calculation, Co is the noise covariance matrix but unknown.
%%%   ..assume it is the identity matrix
Co=zeros(CONST.N);
for i=1:CONST.N; Co(i,i)=1;
end

%do direction finding for each good FOregi
for i4=1:1:length(FOregi)   
    %%% isolate point of interest for this loop
    bp=FOregi(i4,:);   
    %%% make covariance function,
    %loop over all antennas to fill out C
    C=nan.*ones(CONST.N,CONST.N); 
    for ii=1:CONST.N
        eval(['C(ii,ii)=data.a' num2str(ii) num2str(ii)   '(bp(1),bp(2));'])
        for jj=ii+1:CONST.N
            eval(['C(ii,jj)=data.a' num2str(ii) num2str(jj)   '(bp(1),bp(2));'])
            eval(['C(jj,ii)=conj(data.a' num2str(ii) num2str(jj)   '(bp(1),bp(2)));'])
        end
    end
    
    if isempty(find(isnan(C(:))==1))==1  %fix to ensure that data in C is good
        %%%% preprocess the eig function, to be used in music, mle, etc.
        %%%% equation to solve is c = a*s*a*T
        [Q ,D]=eig(C); %Compute eigendecomposition of covariance matrix
        [D,I]=sort(diag(D),1,'descend'); %Find r largest eigenvalues
        Q=Q(:,I); %Sort the eigenvectors to put signal eigenvectors first
        
        
        %%% examine eigenvalues to get statistical determination of the number of signals present
        %%%   a priori, how many solutions there should be... follows Johnson and Dudgen, sec 7.3.5
        M=CONST.N;
        sig_am=nan.*ones(1,M); sig_gm=sig_am;
        MDL=sig_am;  AIC=sig_am;
        L=SpecHead.nsegs.*(SpecHead.PC.SpecOverlapPct/100);
        Ns=CONST.N-1:-1:0;  MDL=nan.*ones(1,CONST.N); AIC=MDL; sig_gm=MDL;  sig_am=MDL;
        for ii=1:length(Ns);
            sig_am(ii)= (1/(CONST.N-Ns(ii))) * sum(D(Ns(ii)+1:CONST.N));  %uses arthmetic mean
            sig_gm(ii)= (prod(D(Ns(ii)+1:CONST.N)))^(1/(CONST.N-Ns(ii)))  ;   %uses geometric mean  (i.e. the nth root of the product of the n values
            MDL(ii)= -L.*(M-Ns(ii)).*log((sig_gm(ii))./sig_am(ii)) +  .5*Ns(ii)*(2*M-Ns(ii)+1)*log(L) ;
            AIC(ii)= -L.*(M-Ns(ii)).*log((sig_gm(ii))./sig_am(ii)) +  Ns(ii)*(2*M-Ns(ii)+1) ;
        end
        
        %%%    find the min value of MDL/AIC as the assumed number of signals present
        Ns_f=nan;
        if strcmp(CONST.which_Ns_meth,'MDL')==1
            j=find(Ns<=CONST.Nmax);    i=find(MDL(j)==min(MDL(j)));    Ns_f=Ns(j(i));
        elseif strcmp(CONST.which_Ns_meth,'AIC')==1
            j=find(Ns<=CONST.Nmax);    i=find(AIC(j)==min(AIC(j)));    Ns_f=Ns(j(i));
        end
        
        if CONST.goplot(2)==1;
            figure(20); clf;
            subplot(2,1,1);
            plot(Ns,(real(AIC)),'bo-'); hg; set(gca,'xdir','reverse')
            plot(Ns,(real(MDL)),'ro-'); hg ; legend('AIC','MDL')
            title([num2str([i4 Ns_f])])            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  proceed with the direction finding method.
        if Ns_f>0;
            
            %%%  evaluate the DOA function to assess how many scatterers there are
            %%% steps to id which MS is the correct one
            %%%
            %%%  send the covariance, the APM, and theta, and the max number of signals, to
            %%%  the DF method in question ...
            %%%  regardless of the methodology, each will return:
            %%%  doa  is the doa fn or the closest thing to it
            %%%  idx   is the emitter locations, given the assumed number of emitters present
            %%%  k    is the max number of iterations used
            %%
           
            switch CONST.which_df
                %%%%%%%%%%%%% MUSIC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 'music'
                    %%
     
                    [doa,idx] = HFR_DP_DFmusic_v1(patt.A,Q,Ns_f,CONST.N,peak_thresh);
                    %make this doa 'look' like the rest
                    % music can return a variable number of peaks b/c of the peak
                    % finder threshold
                    idx=idx{:};
                    [~,i]=sort(doa(idx),'descend'); idx=idx(i);
                    if length(idx)>Ns_f;            idx=idx(1:Ns_f);        end
                    doa=repmat(real(doa),[length(idx),1]);
                    k=nan;
                    
                    %%%%%%%%%%%%% mle  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 'mle'
                    %mle
                    [doa,idx,k] = HFR_DP_DFmle_v1(patt.A,C,Ns_f);
                    %%%%%%%%%%%%% wsf  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 'wsf'
                    [doa,idx,k] = HFR_DP_DFwsf_v1(patt.A,C,Q,D,Ns_f);
                  %  doa=-doa;                   
            end
            
            if CONST.goplot(2)==1;
                subplot(2,1,2);
                %%% plot this result so far %%%
                colors='krbgcym';
                for j=1:length(idx);
                    plot(patt.angles,real(10*log10(doa(j,:)')),[colors(j) '.'],'markersize',4); hg;
                    plot(patt.angles(idx(j)),real(10*log10(doa(j,idx(j)))),[colors(j) 'o'],'markersize',10);
                end
                title(num2str([bp Ns_f])); ylabel('dB')
                pause
            end
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% export the results.
            
            %%%%%% continue to calculate power, and output variables in schmidt notation %%%
            %%% loop over the number of peaks for this solution and export the
            %%% results into one line of RM file for each solution
            
            for iii=1:length(idx)
                %%% form the exporting structure, similar for all parts below
                M=[];
                M.Eigen=abs(D);   % we should give all the eigenvalues every time
                %find mean snr of all antennas
                M.snrant=[mean(diag(C)./noise_floor) ];   %the way to do this if noise_floor is in volts
               % M.snrant=[gain3(bp(1),bp(2))-noise_floor];  %the way to do this if gain3 and noise_floor is already in log10 as dB
                M.isaveout=zeros(1,CONST.N-1);
                M.isaveout(1:Ns_f)=1;
                
                %%% gather peak specific data  %%%%%%%%%
                M.angles=[patt.angles(idx(iii))];
                
                %%% use found peak directions to calculate A and P (in schmidt) or S (in lipa et al)
                %%% export the Schmidt estimated power, as it is different from the measured
                A1=patt.A(:,idx(1:iii));
                %    eval(['P' num2str(ii) '=((A1''*A1)^-1)*A1''*(C- D(length(itheta)+1)*Co)*A1*((A1''*A1)^-1);'])  %back out the signal power
                P=((A1'*A1)^-1)*A1'*(C- D(length(idx)+1)*Co)*A1*((A1'*A1)^-1); %  %back out the signal power
                M.power=P(iii,iii);  %straight up signal voltage?
                %M.power=10*log10(abs(P(iii,iii)))+ (-40 - 5.8); %convert to dB
                M.DOApeak=doa(iii,idx(iii));
                
               
                %%% find the 3dB width of this peak, this can be tricky if the peak is
                %%% small or there are other peaks around, that are larger than this
                %%% one.
                a=3; %3 dB by default
                %%% (1)  adjust cutoff downward if needed
                if range(doa(iii,:))<a   %sometimes the whole thing is quite small, adjust the 3db downwards to be 1/4 of the total range
                    a=max([range(real(doa(iii,:)))./4]);
                end
                l=fliplr(doa(iii,1:idx(iii)));
                r=doa(iii,idx(iii):end);
 
                if strcmp(CONST.which_df,'wsf')==1;  % b/c wsf 'peaks' are local minimums
                %%%  (2) find the next local minimums.
                    li=[find(diff(sign(diff(l(~isnan(l)))))==-2)  length(l)-1];
                    ri=[find(diff(sign(diff(r(~isnan(r)))))==-2)  length(r)-1];
                    %%% (3) look along r and l to see if we can get 3 db up before
                %%% the local min.
                    lii=find(l<l(1)+a);
                    rii=find(r<r(1)+a);
                else
                %%%  (2) find the next local minimums.
                    li=[find(diff(sign(diff(l(~isnan(l)))))==2)  length(l)-1];
                    ri=[find(diff(sign(diff(r(~isnan(r)))))==2)  length(r)-1];
                %%% (3) look along r and l to see if we can get 3 db down before
                     %%% the local min.
                     lii=find(l>l(1)-a);
                     rii=find(r>r(1)-a);             
                end     
                
                %%% (5) pick the smallest,
                M.halfpowwidth(1)=nanmin([length(rii)+length(lii) nanmin([li(1) ri(1)])*2]);
                
                
                M.param=[k Ns_f a];  %output results for the df method
                M.DF_flag=0;
                if length(rii)+length(lii) < nanmin([li(1) ri(1)])*2 % make the flag 1 if using a width from the local minimum
                    M.DF_flag=1;
                    % M.param(2)=1;
                end
                
                if CONST.goplot(2)==1 %& sum(M.isaveout)>2;   M
                    M
                %      figure(4); clf; plot(l); hg; plot(r);
                    pause
                end
                
                % %    M.crit=[nan nan nan];
                % %     %choose the solution, keep this solution and output the required information
                % %     if (7-(M.crit*[1 2 4]'))==0; jj=2;   %use the duel or higher angle solution, follow the DF_flag way of doing this
                % %     else; jj=1; %use the single or lower angle solution.
                % %     end
                
                %%% The flag %%%%
                %%% since the successful number of peak info is elsewhere, total_flag only
                %%% needs to show why a single/duel (for example) was made
                %combine all flags together (each get 2 digits)
                %M.DF_flag=(7-(M.crit*[1 2 4]'));
                %for each iteration of the comparison.
                
                %%% after the method is complete export M in the required format
                %%%% put outgoing data into an array similar to radial metrics, starting
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
                r=[d bp Ns_f iii abs(M.power) M.halfpowwidth M.DOApeak M.snrant M.param M.isaveout M.Eigen' M.DF_flag];
                eval(['R' num2str(iii) '(i4,:)=r;'])
                
                %  Ns_f and iii tell you what solution was chosen and what angle of this
                %  solution is given in the line of data.
                %%%%%%%% end export data section %%%%%%%%%
                               
            end % for iii length(idx) loop over all peaks in the single solution
        end  %if  Ns_f is not empty
    end   %if C has no nans
    
end  % for i4  FOregi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%combine all Rs, while getting rid of empty lines (marked by nans)
for ii=1:CONST.Nmax
    eval(['i' num2str(ii) '=find(isnan(R' num2str(ii) '(:,1))==0);'])
end
if CONST.Nmax==4
    R=[R1(i1,:); R2(i2,:); R3(i3,:); R4(i4,:)];
elseif CONST.Nmax==5
    R=[R1(i1,:); R2(i2,:); R3(i3,:); R4(i4,:); R5(i5,:)];
end
clear R1 R2 R3 R4 R5


%%%%% add this script to the processing steps
HEAD.ProcessingSteps{end+1}=mfilename;

return

