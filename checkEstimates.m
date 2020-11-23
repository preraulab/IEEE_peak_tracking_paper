function measures = checkEstimates(ss_est,ss_mdl,y_true,x_true,alpha_true,ss_gen,K_lim,T_lim)
% checkEstimates takes EKF/IEKF estimates of a multi-peak state-space model
% applied to simulated data and compares them to the true values to 
% evaluate the filter performance based on a set of statistical measures:
%
% It uses several internal functions given below:
%   icTest, Lmatr, Dmatr, combo2PickParams, combo2PickPeaks, getComboNum
%
% INPUTS:
%   ss_est     -- estimate struct containing filter estimate variables.
%                 Required.
%   ss_mdl     -- StateSpaceMultiPeak object used to obtain filter
%                 estimates. Required.
%   y_true     -- matrix of true observations (dim_y x NT). Required.
%   x_true     -- matrix of true states (dim_x x NT+1). Required.
%   alpha_true -- indicator matrix of true On/Off-combinations 
%                 (num_combos x NT+1). Required.
%   ss_gen     -- StateSpaceMultiPeak object used to generate the true
%                 observations. Required.
%   K_lim      -- Positive integer sets number of frequency bins over which 
%                 correlations are evaluated in icTest, boxQ, and moranI. 
%                 Default 3.
%   T_lim      -- Positive integer sets number of time steps over which 
%                 correlations are evaluated in boxQ and moranI. Default 4. 
%
% OUTPUTS:
%   measures -- a vector containing statistical values (15 x 1)
%
% Created by Patrick Stokes
%

%*************************
% Handle variable inputs *
%*************************
if nargin < 8
   T_lim = [];
end
if nargin < 7
   K_lim = [];
end
if nargin < 6
   error('Insufficient inputs provided');
end

%******************************
% Set defaults / Check inputs *
%******************************
if isempty(T_lim)
    T_lim = 4; 
end
if isempty(K_lim)
    K_lim = 3;
end
if isempty(ss_gen)
    error('Input ss_gen was empty. It is required');
end
if isempty(alpha_true)
    error('Input alpha_true was empty. It is required');    
end
if isempty(x_true)
    error('Input x_true was empty. It is required');
end
if isempty(y_true)
    error('Input y_true was empty. It is required');
end
if isempty(ss_mdl)
    error('Input ss_mdl was empty. It is required');
end
if isempty(ss_est)
    error('Input ss_est was empty. It is required');
end

%*********************
% Storage for output *
%*********************
measures = zeros(15,1);

%**********************************************
% Compute errors and other needed information *
%**********************************************
NT = size(y_true,2);
dim_x = ss_mdl.dimX;
idx_t = 2:(NT+1);

% Filter observation error
resid_y = y_true - ss_est.yf_hat(:,idx_t);

% Filter state error 
resid_x = x_true(:,idx_t) - ss_est.xf_hat(:,idx_t);

% Indicator matrix (dim_x x NT) of when parameters are on
pick_on_var_true = logical(combo2PickParams(alpha_true(:,idx_t),ss_gen));

% Indicator matrix (dim_x x NT) of when parameters are on
pick_on_peaks = combo2PickPeaks(ss_est.alpha(:,idx_t),ss_mdl);

% Indicator matrix (num_peaks x NT) of when peaks are on
pick_on_peaks_true = combo2PickPeaks(alpha_true(:,idx_t),ss_gen);

% Vector (1 x NT) of indices of estimated On/Off-combinations
combo_num = getComboNum(ss_est.alpha(:,idx_t));

% Vector (1 x NT) of indices of true On/Off-combinations
combo_num_true = getComboNum(alpha_true(:,idx_t));


%*****************************************
% Compute statistical measures and tests *
%*****************************************

% mean of filter observation error
measures(1) = mean(resid_y(:));

% MSE of filtered spectrogram estimate
measures(2) = mean(resid_y(:).^2);

% median abs observation error
measures(3) = median(abs(resid_y(:)));

% max abs filter observation error
measures(4) = max(abs(resid_y(:)));

% Statistic of instantaneous correlations 
[ic_stat,~,~,ic_test] = icTest(resid_y,[],0,K_lim);
measures(5) = ic_stat;
measures(6) = ic_test;

% Box Q whiteness statistic
[b_stat,~,~,b_test] = boxQ(resid_y,ceil((2*T_lim+1)/2),[],2,0,K_lim);
measures(7) = b_stat;
measures(8) = b_test;

% Moran I spatial autocorrelation statistic
[i_stat,~,~,~,~,i_test] = moranI(resid_y,K_lim,T_lim,[]);
measures(9) = i_stat;
measures(10) = i_test;

% Sum state square filter errors (over times when peaks are truly on) and
% determine times when state filter estimate lies within the confidence interval
sse = 0;
xf_hi = zeros(dim_x,NT); 
xf_lo = zeros(dim_x,NT); 
covered = zeros(dim_x,NT); 
for ii = 1:dim_x
    sse = sse + sum(resid_x(ii,pick_on_var_true(ii,:)).^2);
    xf_hi(ii,:) = ss_est.xf_hat(ii,idx_t)+2*sqrt(squeeze(ss_est.Pf(ii,ii,idx_t))');
    xf_lo(ii,:) = ss_est.xf_hat(ii,idx_t)-2*sqrt(squeeze(ss_est.Pf(ii,ii,idx_t))');
    covered(ii,:) = (x_true(ii,idx_t)<=xf_hi(ii,:)) & (x_true(ii,idx_t)>=xf_lo(ii,:));    
end

% MSE of filtered state estimate (true peak-on times)
measures(11) = sse/sum(sum(pick_on_var_true));

% State coverage probability (all times)
measures(12) = sum(sum(covered))/(NT*dim_x);

% Cohen's kappa of estimated combo
po = sum(combo_num_true==combo_num)/NT;
pe = sum(sum(ss_est.alpha(:,idx_t),2).*sum(alpha_true(:,idx_t),2))/NT^2;
measures(13) = (po-pe)/(1-pe);

% Mean number of correct on/off-peaks
correct_peaks = pick_on_peaks_true==pick_on_peaks;
measures(14) = mean(sum(correct_peaks,1));

% Probability of correct peak-on determination 
measures(15) = sum(sum(correct_peaks))/(size(correct_peaks,1)*size(correct_peaks,2));

end

function [lam_W, pVal, cVal, lam_test] = icTest(u,alpha,f_ver,K_lim)
% icTest tests the instantaneous cross-correlations against the null of 
% being equal to zero using a Wald statistic, following the Lutkepohl
% test of "instantaneous causality". Two versions are provided. The first 
% version uses the matrix computations from Lutkepohl. The second version 
% is a simpler computation of the sum of square correlations. This 
% version can be restrict to limit the spatial extent of the correlations 
% evaluated.
%
% NOTE: The Lutkepohl version is not viable for large dimensions due
%       a Kronecker product in the computation.
%
% See: 
% Lutkepohl - New Introduction to Multiple Time Series Analysis - 
%             Testing for Instantaneous Causality p. 104-105 - eqn. 3.6.12
% 
% INPUTS:
%   u        -- matrix (K x T) of residuals. Required.
%   alpha    -- significance level. Default is 0.05.
%   f_ver    -- binary flag indicating version to compute.
%               1 - Lutkepohl matrix computation.   
%               0 - simple sum of square correlations. (Default)
%   K_lim    -- integer number limiting the number of spatial bins 
%               across which the correlations are evaluated. Not 
%               applicable to matrix computation version. 
%               Default is K-1.
%
% OUTPUTS:
%   lam_W    -- statistic value
%   pVal     -- chi-square p-value
%   cVal     -- chi-square critical value
%   lam_test -- chi-square test result
%
%

%*************************
% Handle variable inputs *
%*************************
if nargin < 4
    K_lim = [];
end
if nargin < 3
    f_ver = [];
end
if nargin < 2
    alpha = [];
end

%***************
% Set defaults *
%***************
if isempty(alpha)
    alpha = 0.05;
end
if isempty(f_ver)
    f_ver = 0;
end

[K,T] = size(u);
if isempty(K_lim)
    K_lim = K-1;
end

%*************************
% Compute Wald statistic *
%*************************
u = u - repmat(mean(u,2),[1 T]);
if f_ver==1
    Sig_u = cov(u');
    L = Lmatr(K);
    sig = L*reshape(Sig_u,[K^2 1]);
    D = Dmatr(K);
    C = zeros(K*(K-1)/2,K*(K+1)/2);
    cnt_c = 0;
    cnt_r = 0;
    for ii = 1:K
        C((cnt_r+1):(cnt_r+K-ii),(cnt_c+2):(cnt_c+K-(ii-1))) = eye(K-ii);
        cnt_c = cnt_c+K-(ii-1);
        cnt_r = cnt_r+K-ii;
    end
    term = 2*C*pinv(D)*kron(Sig_u,Sig_u)*pinv(D)'*C';
    lam_W = T*sig'*C'*(term\C*sig);
else
    C = zeros(K,K);
    for tt = 1:T
        C = C + u(:,tt)*u(:,tt)';
    end
    C = C/T;
    C = C./(repmat(sqrt(diag(C)),[1 K]).*repmat(sqrt(diag(C))',[K 1]));
    if K_lim < (K-1)
        C = tril(triu(C,-K_lim),K_lim);
    end
    C = tril(C,-1);
    C = C.*C;
    lam_W = T*sum(C(:));
end

%**************************
% Perform chi square test *
%**************************
dof = K*(K-1)/2;
if K_lim < (K-1)
    dof = dof - (K-K_lim)*(K-K_lim-1)/2;
end
pVal =  chi2cdf(lam_W,dof,'upper');
cVal = chi2inv(1-alpha,dof);
lam_test = lam_W > cVal;

end
function L = Lmatr(m)
% Forms the elimination matrix L that converts vec (the column-wise 
% stacking of (mxm) matrix elements into a vector) into vech (the column-wise 
% stacking of lower triangular elements). For symmetric matrices,
% the inverse is given by the duplication matrix. Used in icTest.
L = zeros(m*(m+1)/2,m^2);
Im = eye(m);
cnt = 0;
for ii = 1:m
    L((cnt+1):(cnt+m-(ii-1)),(1+(ii-1)*m):ii*m) = Im(ii:end,:);
    cnt = cnt + m-(ii-1);
end
end
function D = Dmatr(m)
% Forms the duplication matrix D that converts vech (the column-wise 
% stacking of lower triangular elements of a (mxm) symmetric matrix) 
% into vec (the column-wise stacking of all matrix elements into a vector). 
% The inverse is given by the elimination matrix. Used in icTest.
D = zeros(m^2,m*(m+1)/2);
cnt = 0;
for ii = 1:m
    cnt2 = 0;
    for jj = 1:(ii-1)
        idx = cnt2 + ii;
        D(jj+(ii-1)*m,idx) = 1;
        cnt2 = cnt2 + m-(jj-1)-1;
    end
    D((ii+(ii-1)*m):ii*m,(cnt+1):(cnt+m-(ii-1))) = eye(m-(ii-1));
    cnt = cnt + m-(ii-1);
end
end

function pick_on_params = combo2PickParams(pick_combo,m)
% Converts the indicator matrix of combos into an indicator matrix of 
% parameters. 
% INPUT: pick_combo -- indicator matrix (num_combos x NT) of
%                      On-combination at each time.
% OUTPUT: pick_on_params -- indicator matrix (num_params x NT) of 
%                           On-peaks at each time.
%
N = size(pick_combo,2);
pick_on_params = zeros(m.totalPeakParams,N);
for nn = 1:N
    tmp_combo = pick_combo(:,nn)>0;
    % Find which peaks are on
    if sum(tmp_combo)>1
        % more than one combo has positive probability
        % take the union of all such peaks
        pick_peaks = sum(m.peakCombos(tmp_combo,:),1)>0;
    else
        pick_peaks = m.peakCombos(tmp_combo,:);
    end
    % Find parameters of On peaks
    pick_on_params(:,nn) = m.getPeakIdxs(logical(pick_peaks));
end
end


function pick_on_peaks = combo2PickPeaks(pick_combo,m)
% Converts the indicator matrix of combos into an indicator matrix of peaks. 
% INPUT: pick_combo -- indicator matrix (num_combos x NT) of
%                      On-combination at each time.
% OUTPUT: pick_on_peaks -- indicator matrix (num_peaks x NT) of 
%                          On-peaks at each time.
%
N = size(pick_combo,2);
pick_on_peaks = zeros(m.numPeaks,N);
for nn = 1:N
    tmp_combo = pick_combo(:,nn)>0;
    if sum(tmp_combo)>1
        % more than one combo has positive probability
        % take the union of all such peaks
        pick_on_peaks(:,nn) = sum(m.peakCombos(tmp_combo,:),1)>0;
    else
        pick_on_peaks(:,nn) = m.peakCombos(tmp_combo,:);
    end
end
end

function combo_num = getComboNum(pick_combo)
% Converts the indicator matrix of combos into a vector 
% of combo indices.
% INPUT: pick_combo -- indicator matrix (num_combos x NT) of
%                      On-combination at each time.
% OUTPUT: combo_num -- vector (1 x NT) of On-combination index at each
%                      time.
%
combo_num = zeros(1,size(pick_combo,2));
for ii = 1:length(combo_num)
    combo_num(ii) = find(pick_combo(:,ii));
end

end

