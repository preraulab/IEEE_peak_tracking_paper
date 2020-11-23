function [Y, Xbar, alpha, is_on, V, comps] = simulatePseudoDeterm(ss,stimes,rand_seed)
% simulatePseudoDeterm generates a simulation of a multi-peak state-space
% model except with parameters following randomly determined patterns:
% cosines (40%), saw-tooths (20%), steps (30%), and exponential cusps (10%).
%
% INPUTS:
%   ss        -- StateSpaceMultiPeak object containing the peaks models,
%                parameter bounds, observation noise covariance, etc.
%                Required.
%   stimes    -- vector (1 x N) of times at which the simulation is
%                evaluated. Default is 0.005*(1:500).
%   rand_seed -- random seed to recreate simulation. Default is empty [].
%
% OUTPUTS:
%   Y     -- simulated spectrogram matrix (dim_Y x N)
%   Xbar  -- matrix (dim_X x N) of simulated unbounded states
%   alpha -- indicator matrix (num_combos x N) of On/Off-combination
%   is_on -- indicator matrix (num_peaks x N) of On/Off-peaks
%   V     -- simulated observation noise matrix (dim_Y x N)
%   comps -- array (num_peaks x dim_Y x N) of simulated individual peaks
%
% Created by Patrick Stokes
% 

%*************************
% Handle variable inputs *
%*************************
if nargin < 3
    rand_seed = [];
end
if nargin < 2
    stimes = [];
end
if nargin < 1
    error('StateSpaceMultiPeak object ss must be provided.');
end

%***************
% Set defaults *
%***************
if isempty(stimes)
    dt = .005;
    N = 500;
    stimes = dt*(1:N); 
else
    N = length(stimes);
    dt = stimes(2)-stimes(1);
end
if ~isempty(rand_seed)
    rng(rand_seed);
end
if isempty(ss)
    error('StateSpaceMultiPeak object ss must be provided.');
end

%***************************
% Deterministic parameters *
%***************************
% Random pattern type for each parameter
sim_funcs = randi(10,ss.totalPeakParams,1);
% Random frequency of patterns (if used)
sim_freqs = rand(ss.totalPeakParams,1)*(7/max(stimes))+1/max(stimes);
% Random phase of pattern (if used)
sim_phases = rand(ss.totalPeakParams,1)*2*pi;

%************************************
% Form X (bounded parameter values) *
%************************************
X = zeros(ss.totalPeakParams,N);
cnt_param = 1;
for ii = 1:ss.numPeaks
    tmp_bnds = {ss.peakModels(ii).modelParams.linkFunParams};
    for jj = 1:ss.peakModels(ii).numPeakParams
        tmp_bnd = tmp_bnds{jj};
        tmp_rng = tmp_bnd(2)-tmp_bnd(1);
        tmp_lo = tmp_bnd(1)+0.05*tmp_rng;
        tmp_hi = tmp_bnd(2)-0.05*tmp_rng;
        tmp_rng = 0.9*tmp_rng;
        tmp_mid = tmp_lo+tmp_rng/2;
        if sim_funcs(cnt_param) <= 4 % cosine
                tmp = (tmp_rng/2)*cos(2*pi*sim_freqs(cnt_param)*stimes + sim_phases(cnt_param))+tmp_mid;
        elseif sim_funcs(cnt_param) <= 6 % saw-tooth
                tmp = tmp_rng*mod(stimes,1/sim_freqs(cnt_param))*sim_freqs(cnt_param)+tmp_lo;
        elseif sim_funcs(cnt_param) <= 9 % step
                tmp_step = cumsum(rand(5,1));
                tmp_step = max(stimes)*tmp_step/max(tmp_step);
                tmp_amp = tmp_rng*rand(5,1);
                tmp = (stimes<=tmp_step(1))*(tmp_amp(1)+tmp_lo);
                for kk = 2:5
                    tmp(stimes>tmp_step(kk-1) & stimes<=tmp_step(kk)) = tmp_amp(kk)+tmp_lo;
                end
        else % exp
                tmp = [tmp_rng*exp(stimes(1:ceil(N/2))-stimes(ceil(N/2))) tmp_rng*exp(stimes(ceil(N/2))-stimes((ceil(N/2)+1):N))] + tmp_lo;
        end
        X(cnt_param,:) = tmp;
        cnt_param = cnt_param + 1;        
    end
end

%*****************************************************
% Convert X to Xbar (invert link to unbounded state) *
%*****************************************************
Xbar = X;
cnt_param = 1;
for ii = 1:ss.numPeaks
    for jj = 1:ss.peakModels(ii).numPeakParams
        Xbar(cnt_param,:) = ss.peakModels(ii).modelParams(jj).linkInvFun(X(cnt_param,:));
        cnt_param = cnt_param + 1;
    end
end

%**************************
% Determine peak-on times *
%**************************
% On-off times of each peak are determined separately by generating
% intervals as exponential RVs.
is_on = zeros(ss.numPeaks,N);
is_on(:,1) = 1*(rand(ss.numPeaks,1)>0.5);
E_time_onoff = max(stimes)/6;
for ii = 1:ss.numPeaks
    idx_curr = 1;
    t_curr = stimes(idx_curr);
    while t_curr < max(stimes)
        dt_switch = exprnd(E_time_onoff);
        idx_switch = find(stimes<(t_curr+dt_switch),1,'last');
        is_on(ii,idx_curr:idx_switch) = is_on(ii,idx_curr);
        idx_curr = idx_switch+1;
        if idx_curr <= length(stimes)
            is_on(ii,idx_curr) = abs(is_on(ii,idx_curr-1)-1);
            t_curr = stimes(idx_curr);
        else
            t_curr = max(stimes);
        end
    end    
end
is_on(~ss.isDynamic,:) = 1;
% Determine combos from which peaks are on
alpha = zeros(ss.numCombos,N);
for tt = 1:N
    pick_combo = sum(repmat(is_on(:,tt)',[ss.numCombos 1]) == ss.peakCombos,2) == ss.numPeaks;
    alpha(pick_combo,tt) = 1;
end

%********************
% Form observations *
%********************
Y = zeros(ss.dimY,N);
comps = zeros(ss.numPeaks,ss.dimY,N);
V = randn(N,ss.dimY)*chol(ss.R(1));
V = V';
for tt = 1:N
    [~, comps(:,:,tt)] = ss.getH(Xbar(:,tt),ss.omega,is_on(:,tt));
    Y(:,tt) = comps(:,:,tt)'*is_on(:,tt) + V(:,tt);
end
