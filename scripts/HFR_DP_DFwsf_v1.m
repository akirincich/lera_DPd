function [doa,idx,k] = HFR_DP_DFwsf_v1(A,R,U,D,q)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [doa,idx] = HFR_DP_DFwsf_v1(A,R,U,D,q)

% WEIGHTED SUBSPACE FITTING -  Estimate DOA with WSF algorithm
%
% INPUTS
% A    - Array Manifold (aka array matrix at all theta (complex))
% R    - Covariance matrix
% U    - the eigevectors of R, ordered
% D    - the eigenvalues of R, ordered
% q    - max number of signal sources (must be max of N-1)
%
% *If a radar has M elements and theta has d bearings, A is Mxd
%
%
% OUTPUTS
% doa    - the direction of arrival function for q assumed signals
% idx    - Cell array with indicies of the q signal source solutions
%          eg idx{1} has the index of the single source soln, idx{2} is
%          dual, etc
%
% REFERENCE
% Krim and Viberg, 1996, "Two Decades of Array Signal Processing Research: 
%   the Parametric Approach", IEEE
% Viberg, Ottersten and Kailath, 1991, "Detection and Estimation in Sensor 
%   Arrays Using Weighted Subspace Fitting", IEEE  
%

% adapted from
% Copyright (C) 2017 Brian Emery
%
% Version 13-Jun-2017 
%
%v1  adapted from emery to also give the doa back
% Anthony Kirincich 
% WHOI-PO
% akirincich@whoi.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%% test 
%A=patt.A; R=C;  U=Q; q=3;

%%
% INITIALIZATION - use MLE-AP for first pass 
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

%%%%%%%%%%%%%%%%%%%% INIT WSF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use sorted eigen decomposition of covariance matrix input

% Noise and Signal sub spaces for n brg soln
Us = U(:,1:q);
% Matrix of Signal eigenvalues
Ds = diag(D(1:q)); %D is presorted.

% estimate noise variance from noise eigen values % < --- CHECK
% does this need to be squared?
s = mean(D(q+1:end));
% Compute the Weight matrix
W = (( Ds - s.*eye(q) )^2) /(Ds); % /(Ds) is eqivalent to *inv(Ds)

%%%%%%%%%%%%% MAIN LOOP for WSF  %%%%%%%%%%%%%%%%%%%%%%
d = ones(size(q)); 
k = 1;
thi_k = thi; % thi is current iteration, thi_k is the next

doa = nan.*ones(q,n); %save the trace result for the final run.
I = eye(size(A,1));    % form an identity matrix with size Nmax

while any(d > eps) 
    for jj = 1:q
        % index of thi to use in prior projection
        ith =  setdiff(1:q,jj);          
        % loop over the rest of thi, for this iteration
% thi_k(i) = arg_min(A,Us,W,thi(x));
% function thi = arg_min(A,Us,W,thx)   % Simple version of MLE-AP ... for use with WSF.M

% inti storage for this iteration
tr = NaN(n,1);
for ii = setdiff(1:n,thi(ith))
    % Projection matrix (function of theta) PI_A perp
    PI = I - (A(:,[thi(ith) ii])*pinv(A(:,[thi(ith) ii])));   
    % compute trace
    tr(ii) = trace( PI*Us*W*Us' );  
end   % loop over all array manifold values (directions) in the
[~,thi_k(jj)] = min(tr);

doa(jj,:)=real(tr);   %save the trace in case this is the final time
    
        % update this to reflect the change in thi 
        thi = thi_k;
    end   % loop over jj=1:q 

         % check for convergence, with this jj change in thi_k
        d = sqrt(( thi - thi_k ).^2);
        k = k+1;  % disp(num2str(k))   

    % Put a break in here?
    if k == 50,
        disp('WSF: 100 iterations reached, terminated')
        break
    end
end   %while

idx = thi(:);

return


