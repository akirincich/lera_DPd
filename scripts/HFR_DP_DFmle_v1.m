function  [doa,idx,k] = HFR_DP_DFmle_v1(A,R,q);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function    [doa,idx,k] = HFR_DP_DFmle_v1(A,R,q);
%
% input variables are:
%   A                    Array Manifold (aka array matrix at all theta
%                        (complex)),  If a radar has M elements and
%                           theta has d bearings, A is Mxd
%   R  (or C)           the observed covariance matrix
%   q                   # of assumed incoming signals
%
% output variables are:
%   doa        the DOA functions for each of the q solutions
%   idx        the indices of the DOA funtion that are defined peaks
%   k          the max number of iterations made to define the signal
%               locations
%
%
%  Note that this code was adapted directly from B.Emery's mle_ap.m 
%   but condensed to skip test functionality as well as return trace
%   results.
%
% 11/2017
%Anthony Kirincich
%akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% %%% test inputs, have A,C,ii?
% A=patt.A; q=3;  R=C;

% INITIALIZATION,   Run through the projection and trace for each emitter
% choice separately... using eqn 17 and 18, etc. from
thi = [];
n = size(A,2);       % Get number of bearings to loop over for the maximization

for ii = 1:q
    % get the initial indicies of the q values for theta_i,
    B = A(:,thi);   %even if thi is empty, this forms correctly
    tr = NaN(n,1);      % init storage of trace output, for this ii run
    
    % loop over index of possible, excluding known indices, thi
    for jj = setdiff(1:n,thi)
        a=[B A(:,jj)];   % get matrix of Array manifold for known indices
        %and the indix in question    %P = proj([B A(:,i)]);
        % COMPUTE PROJECTION MATRIX, Note: ' is the Hermitian transpose, which is intended
        P = a*( (a'*a)\a' );
        tr(jj) = trace(P*R);   %store the trace, the sum of the eigenvalues
    end   %for jj loop over values of the array manifold.
    [~,ix] = max(tr);   %find the simple peak.
    thi = [thi ix];        % add them in to make the augmented matrix
end   %for ii loop over q possible emitters
%%

% MAIN LOOP, do this again, using only the other emitters in the projection
% to refine the estimate of the index of the emitter in quesion
d = ones(size(q)); 
k = 1;    %set up for continuing loop while thi are still changing.
thi_k = thi; % thi is current iteration, thi_k is the next

doa = nan.*ones(q,n); %save the trace result for the final run.

while any(d > eps)    %d is non-zero
    for jj = 1:q        
        % index of thi to use in prior projection
        ith =  setdiff(1:q,jj);
        % compute the projection for the non-jj sources,  %Pb = proj(A(:, thi(ith) ));
        a=A(:, thi(ith) );
        % COMPUTE prior PROJECTION MATRIX, Note: ' is the Hermitian transpose, which is intended
        Pb = a*( (a'*a)\a' );
        % form an identity matrix with size pb
          I = eye(size(Pb));    
        % loop over all values of thi, for this iteration  %thi_k(i) = arg_max_b(A,R,Pb,thi(x));               
        % this will form the augmentend matrix for multi-bearing cases

        tr = NaN(n,1);  %storage of the current trace values.
        for ii = setdiff(1:n,thi(ith))   %for all values but those indices included in a and Pb
            % % update the projection by forming an identidy matrix %b = proj_update(A(:,i),Pb);
            %  I = eye(size(Pb));    
            % update the projection by forming b %b = proj_update(A(:,i),Pb);
            % 'a(theta) sub A(THETA)'
            Cb = (I - Pb)*A(:,ii);
            b = Cb/norm(Cb); % eqn 22.b
            tr(ii) = b'*R*b;    %and the updated trace...
        end   % loop over all array manifold values (directions) in the 
        [~,thi_k(jj)] = max(tr);    %find the simple max
        
        doa(jj,:)=real(tr);   %save the trace in case this is the final time
        
        % check for convergence, with this jj change in thi_k
        d = sqrt(( thi - thi_k ).^2);
        % update thi to reflect 
        thi = thi_k;
        k = k+1;  % disp(num2str(k))   
    end   % loop over jj=1:q 
    
    % Put a break in here?
    if k == 100,
        disp('MLE: 100 iterations reached, terminated')
        break
    end 

end    %while

idx = thi(:);  % push to the output


return
