function [ tmatrix ] = makeComboTransitionMatr(peak_combos, p_transition)
% makeComboTransitionMatr determines the combo transition matrix for a
% given set of peak combos from a set of input probabilities. 
% 
% INPUTS
%   peak_combos  -- indicator matrix (num_combos x num_peaks) of which
%                   peaks are on in each combo. Required.
%   p_transition -- Four possibilities:
%                   [p1]       - peaks turn on and off with probability p1
%                   [p1 p2]    - peaks turn on with probability p1
%                                and turn off with probability p2
%                   [p1 p2 p3] - probability of staying in the current
%                                combo is p3 and transitioning is (1-p3).
%                                Away transitions determined by
%                                probabilities of peaks turning on (p1) and
%                                turning off (p2).
%                   otherwise  - uniform transition probabilities between
%                                combos
%                   After these setups, the overall transition 
%                   probabilities are then appropriately normalized. 
%                   
% OUTPUTS:
%   tmatrix -- matrix (num_combos x num_combos) of transition probabilities
%
% Created by Patrick Stokes
%

%*************************
% Handle variable inputs *
%*************************
if nargin < 2
    p_transition = [];
end
if nargin < 1 || isempty(peak_combos)
    error('Matrix of peak combos, peak_combos, required.');
end

%*************************
% Determine input option *
%*************************
num_combos = size(peak_combos,1);
% Initializing with ones handles uniform case
tmatrix = ones(num_combos,num_combos);

if length(p_transition)==1
    p_on = p_transition;
    p_off = p_transition;
    f_unif = false;
elseif length(p_transition)>=2
    p_on = p_transition(1);
    p_off = p_transition(2);
    f_unif = false;
else
    f_unif = true;
end
if length(p_transition)==3
    p_stay = p_transition(3);
else
    p_stay = 0;
end

%*************************
% Form transition matrix *
%*************************
% Set "relative probabilities" based p_on and p_off
if ~f_unif
    for ii=1:size(peak_combos,1)
        for jj=1:size(peak_combos,1)
            % tmatrix(ii,jj)=p_transition^sum(~(peak_combos(ii,:) == peak_combos(jj,:)))*(1-p_transition)^sum((peak_combos(ii,:) == peak_combos(jj,:)));
            for kk = 1:size(peak_combos,2)
                if peak_combos(jj,kk)==0
                    tmatrix(ii,jj) = tmatrix(ii,jj)*p_on^(~(peak_combos(ii,kk) == peak_combos(jj,kk)))*(1-p_on)^((peak_combos(ii,kk) == peak_combos(jj,kk)));
                elseif peak_combos(jj,kk)==1
                    tmatrix(ii,jj) = tmatrix(ii,jj)*p_off^(~(peak_combos(ii,kk) == peak_combos(jj,kk)))*(1-p_off)^((peak_combos(ii,kk) == peak_combos(jj,kk)));
                else
                    
                end
            end
        end
    end
end
% Case of uniform transition handled b/c tmatrix was initialized to ones

% Normalize each column
for jj=1:num_combos
    if p_stay~=0
        % Set diagonal element of column to p_stay  
        tmatrix(jj,jj) = p_stay;
        % Normalize off-diagonal elements to sum to (1-p_stay)
        pick_offdiag = setdiff(1:num_combos,jj);
        tmatrix(pick_offdiag,jj) = tmatrix(pick_offdiag,jj).*(1-p_stay)./sum(tmatrix(pick_offdiag,jj));
    else
        tmatrix(:,jj) = tmatrix(:,jj)./sum(tmatrix(:,jj));
    end
end

end
