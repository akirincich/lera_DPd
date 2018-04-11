function    [music_spectrum,idx] = music_v1(A,Q,q,N,peak_thresh);
    
%%%%%%%%%%%%%%%%%
%function    [music_spectrum,idx] = music_v1(A,Q,j,N,peak_thresh);
%
% input variables are:
%   A
%   Q                   the ordered eigenvecter matrix
%   q                   # of assumed incoming signals
%   N                   # of antennas
%   peak_thresh         the log-space threshold used to determine a peak in
%                            the DOA fn
%
% output variables are:
%   music_spectrum       the DOA function for j incoming input signals 
%   idx                  the indices of the DOA funtion that are defined peaks
%
%
% 11/2017
%Anthony Kirincich
%akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        Qn=Q(:,q+1:N); %Get the noise eigenvectors       
        music_spectrum=nan.*real(A(1,:));
        for i=1:length(A)
            music_spectrum(i)=(  (A(:,i)'*A(:,i))/(A(:,i)'*(Qn*Qn')*A(:,i)) );           
%           %%% Compute the DOA function (Schmidt 1986, equation 6)
%           DOA(j,i) = 1/(A(:,i)'*(Qn*Qn')*A(:,i)); %does same above, just different by a factor of ~9 for an ideal ant
        end

    %%% get the peaks of all MS, ignore all peaks on the boundary or at the edge of 360 pattern
       [pks,dzdt]=pksfinder(10*log10(real(music_spectrum)),peak_thresh);
        %focus on the upper, isolated peaks
        idx={find(diff(dzdt)==-2)+1}  ;
        
       
return

