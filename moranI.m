function [I,E_I,var_I,pVal,cVal,Itest] = moranI(u,K_lim,T_lim,alpha)
% moranI computes the Moran's I statistic. The statistic is modified from 
% its full form by allowing a restriction of the two-dimensional extent 
% over which the statistic is evaluated. 
%
% See https://en.wikipedia.org/wiki/Moran%27s_I
% 
% INPUTS:
%   u        -- matrix (K x T) of residuals. Required.
%   K_lim    -- integer number limiting the number of spatial bins 
%               across which the correlations are evaluated. 
%               Default is K-1.
%   T_lim    -- integer number limiting the number of time lags 
%               over which the correlations are evaluated. Default is T-1.
%   alpha    -- significance level. Default is 0.05.
%
% OUTPUTS:
%   I     -- I statistic value
%   E_I   -- expected value under the null hypothesis  
%   var_I -- variance under the null hypothesis
%   pVal  -- z p-value
%   cVal  -- z critical value
%   Itest -- z-test result
%
%
% Created: Patrick Stokes
%

%*************************
% Handle variable inputs *
%*************************
if nargin < 4
    alpha = [];
end
if nargin < 3
    T_lim = [];
end
if nargin < 2
    K_lim = [];
end

%***************
% Set defaults *
%***************
if isempty(alpha)
    alpha = 0.05;
end

[K, T] = size(u);
if isempty(K_lim)
    K_lim = K-1;
end
if isempty(T_lim)
    T_lim = T-1;
end

%**********************
% Compute I statistic *
%**********************
u = u - mean(u(:));
% C = xcorr2(u);
% term = sum(sum(C((K-K_lim):(K+K_lim),(T-T_lim):(T+T_lim)),2),1) - C(K,T);
% C0 = C(K,T);
N = K*T;
Wmatr = zeros(K*T);
term = 0;
C0 = sum(u(:).^2);
for ii = 1:K
    for jj = 1:T
        for kk = max(1,(ii-K_lim)):min((ii+K_lim),K)
            for ll = max(1,(jj-T_lim)):min((jj+T_lim),T)
                if ii~=kk || jj~=ll
                    Wmatr(ii+(jj-1)*K, kk+(ll-1)*K) = 1;
                    term = term + u(ii,jj)*u(kk,ll);
                end
            end
        end
    end
end
% Wmatr = Wmatr./repmat(sum(Wmatr,2),[1 K*T]);
term = sum(sum(Wmatr.*(u(:)*u(:)')));
W = sum(Wmatr(:));
I = N/W * term/C0;

%*******************************************
% Compute null expected value and variance *
%*******************************************
E_I = -1/(N-1);
S1 = (1/2)*sum(sum((Wmatr+Wmatr').^2)); % (1/2) * sum_i * sum_j (w_ij + w_ji)^2
S2 = sum((sum(Wmatr,2) + sum(Wmatr,1)').^2,1); % sum_i (sum_j w_ij + sum_j w_ji)^2
S3 = (sum(u(:).^4)/N)/((sum(u(:).^2)/N)^2);
S4 = (N^2-3*N+3)*S1 - N*S2 + 3*W^2;
S5 = (N^2-N)*S1 - 2*N*S2 + 6*W^2;
var_I = (N*S4-S3*S5)/((N-1)*(N-2)*(N-3)*W^2) - E_I^2;

%*****************
% Perform z-test *
%*****************
z = (I-E_I)/sqrt(var_I);
pVal = 2*normcdf(-abs(z));
cVal = norminv(1-alpha);
if z > cVal
    Itest = 1;
elseif z < -cVal
    Itest = -1;
else 
    Itest = 0;
end

end
