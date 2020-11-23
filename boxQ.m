function [Q,pVal,cVal,Qtest] = boxQ(u,h,alpha,f_ver,f_comput,K_lim)
% boxQ computes the vector form of the Box Q (also called Ljung-Box or Box-Pierce) 
% statistic. The statistic is modified from its full form by allowing a restriction
% on the spatial extent of the statistic. Several proposed normalization 
% options are available. 
%
% See: 
% Lutkepohl - New Introduction to Multiple Time Series Analysis - 
%             Portmanteau Tests p. 169-171 - eqn. 4.4.23 and 4.4.24
% https://en.wikipedia.org/wiki/Ljung%E2%80%93Box_test
% 
% INPUTS:
%   u        -- matrix (K x T) of residuals. Required.
%   h        -- integer number of time lags evaluated. Default is 10.
%   alpha    -- significance level. Default is 0.05.
%   f_ver    -- integer flag indicating type of normalization used.
%               1 - Lutkepohl version, T^2/(T-i) (Default)
%               2 - Ljung-Box version, T*(T+2)/(T-i)
%               o.w. - Box-Pierce version, T.
%   f_comput -- binary flag idication whether to compute the statistic  
%               0 - element-wise via for loops
%               1 - using matrix operations (Default).
%   K_lim    -- integer number limiting the number of spatial bins 
%               across which the correlations are evaluated. 
%               Default is K-1.
%
% OUTPUTS:
%   Q     -- Q statistic value
%   pVal  -- chi-square p-value
%   cVal  -- chi-square critical value
%   Qtest -- chi-square test result
%
%
% Created: Patrick Stokes
%

%*************************
% Handle variable inputs *
%*************************
if nargin < 6
    K_lim = [];
end
if nargin < 5 
    f_comput = [];
end
if nargin < 4
    f_ver = [];
end
if nargin < 3
    alpha = [];
end
if nargin < 2
    h = [];
end

%***************
% Set defaults *
%***************
if isempty(h)
    h = 10;
end
if isempty(alpha)
    alpha = 0.05;
end
if isempty(f_ver)
    f_ver = 1;
end
if isempty(f_comput)
    f_comput = 1;
end

[K, T] = size(u);
if isempty(K_lim)
    K_lim = K-1;
end

%**********************************************
% Compute lagged residual covariance matrices *
%**********************************************
u = u - repmat(mean(u,2),[1 T]);
C = zeros(K,K,h+1);
for ii = 0:h
    for tt = (ii+1):T
        C(:,:,ii+1) = C(:,:,ii+1) + u(:,tt)*u(:,tt-ii)';
    end
    % Limits spatial extent to be evaluated 
    if K_lim < (K-1)
        C(:,:,ii+1) = tril(triu(C(:,:,ii+1),-K_lim),K_lim);
    end
end
C = C/T;

%********************
% Compute statistic *
%********************
Q = 0;
% Sum squares of correlations
for ii = 1:h
    term = 0;
    if f_comput == 1
        term = trace(C(:,:,ii+1)'*pinv(C(:,:,1))*C(:,:,ii+1)*pinv(C(:,:,1)));
    else
        for rr = 1:K
            for ss = 1:K
                term = term + (C(rr,ss,ii+1)/sqrt(C(rr,rr,1)*C(ss,ss,1)))^2;
            end
        end
    end
    if (f_ver == 1) || (f_ver == 2)
        term = term/(T-ii);
    end
    Q = Q + term;
end

% Apply normalization
if f_ver == 1 % Lutkepohl version
    Q = Q*T^2;
elseif f_ver == 2 % Ljung-Box version
    Q = Q*T*(T+2);
else % Box-Pierce version
    Q = Q*T;
end

% Check statistic is positive
if Q<0
    warning('oops: Q-statistic negative.');
end

%**************************
% Perform chi square test *
%**************************
if K_lim < (K-1)
    dof = h* ((2*K_lim+1)*(K-K_lim-1)+(K_lim+1)^2);
else
    dof = h*K^2;
end
pVal =  chi2cdf(Q,dof,'upper');
cVal = chi2inv(1-alpha,dof);
Qtest = Q > cVal;

end
