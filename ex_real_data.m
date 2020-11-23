%ex_real_data
% This script is forked from load_data_form_model_v_sleep_ekf_modif.m,
% which is itself based on load_data_form_model_v_sleep_ekf.m 
% It is forked to compute the mt spectrogram (instead of using the
% pre-saved file) and use the new artifact detector.
% It is intended as an option within the test_real_data.m script. 
%
% Modified: 20200624 -- changed to load pre-saved spectrogram
% Created: 20200224 -- Patrick Stokes
%
%

save_fname = []; 

%***************
% Process data *
%***************
load('sop_spect_mt.mat');

% Store original data prior to subselection
spect_trunc = spect;
spect = nanpow2db(spect_trunc);
sfreqs_trunc = sfreqs;
stimes_trunc = stimes;

% Indices for subselecting frequency bins
freq_inds = [find(sfreqs>0,1,'first'):1:find(sfreqs>1.5,1,'first') find(sfreqs>1.5,1,'first')+1:1:find(sfreqs>35,1,'first') find(sfreqs>35,1,'first')+1:5:find(sfreqs>55,1,'first') find(sfreqs>55,1,'first')+1:2:find(sfreqs>65,1,'first') find(sfreqs>65,1,'first')+1:5:length(sfreqs)];

% Indicies for subselecting time points
time_inds = 1:length(stimes);

% Indicator of time points corresponding to artifacts
is_artifact = sum(~isnan(spect),1)==0;

% Subselect data 
sfreqs = sfreqs(freq_inds);
stimes = stimes(time_inds);
spect = spect(freq_inds,time_inds);

% figure; imagesc(stimes,sfreqs_trunc,nanpow2db(spect_trunc)); axis xy; climscale; colormap(jet);

%******************************
% Set up and form peak models *
%******************************
% Types of bounds for each peak amplitude
% 0 = exponential
% 1 = sigmoid
% 2 = abs
f_link_fun_amp = [1 1 1 1 1 1];

% Set the minimum for exponential bounded peak amplitudes  
min_amp = 1;

% Peak models
if f_link_fun_amp(1) == 1
    p_60 = PeakModel('box',{'sigmoid','sigmoid','sigmoid'},{[55 65],[min_amp 50],[.5 10]},6,'60 Hz');
else
    p_60 = PeakModel('box',{'sigmoid','exp','sigmoid'},{[55 65],[1 1],[.5 10]},6,'60 Hz');
end
if f_link_fun_amp(2) == 1
    p_delta = PeakModel('skewgamma',{'sigmoid','sigmoid','sigmoid','sigmoid'},{[1 8],[min_amp 10],[4 50],[.5 1]},[], 'Delta-Theta');
else
    p_delta = PeakModel('skewgamma',{'sigmoid','exp','sigmoid','sigmoid'},{[1 8],[1 1],[4 50],[.5 1]},[], 'Delta-Theta');
end
if f_link_fun_amp(3) == 1
    p_alpha = PeakModel('harmonic', {'sigmoid','sigmoid','sigmoid','sigmoid'},{[8 12],[5 50],[0.5 3],[0 .7]}, 2, 'Alpha+Harm');
else
    p_alpha = PeakModel('harmonic', {'sigmoid','exp','sigmoid','sigmoid'},{[8 12],[1 1],[0.5 3],[0 .7]}, 2, 'Alpha+Harm');
end
if f_link_fun_amp(4) == 1
    p_sigma = PeakModel('gaussian',{'sigmoid','sigmoid','sigmoid'},{[12 16],[min_amp 50],[0.5 3]},[], 'Sigma');
else
    p_sigma = PeakModel('gaussian',{'sigmoid','exp','sigmoid'},{[12 16],[1 1],[0.5 3]},[], 'Sigma');
end
if f_link_fun_amp(5) == 1
    p_slow = PeakModel('skewgamma',{'sigmoid','sigmoid','sigmoid','sigmoid'},{[0 1.5],[min_amp 50],[0.1 5],[.5 1.9]},[], 'Slow');
else
    p_slow = PeakModel('skewgamma',{'sigmoid','exp','sigmoid','sigmoid'},{[0 1.5],[1 1],[0.1 5],[.5 1.9]},[], 'Slow');
end
if f_link_fun_amp(6) == 1
    p_base = PeakModel('expdecay', {'sigmoid','sigmoid','identity'},{[0 50],[0 .5],[]},[],'Baseline');
else
    p_base = PeakModel('expdecay', {'exp','sigmoid','identity'},{[1 1],[0 .5],[]},[],'Baseline');
end
peaks = [p_60 p_delta p_alpha p_sigma p_slow p_base ];

%***************************
% Set up state-space model *
%***************************
% Set baseline and 60Hz to be always on
is_dynamic = [0 1 1 1 1 0]; 

% Peak combos
%            60 d ah sig s bl
peak_combos = [1 0 1 0 1 1; ... % wake
               1 0 1 0 0 1; ... % wake
               1 0 0 0 0 1; ... % n1
               1 0 0 0 1 1; ... % n1
               1 1 0 1 1 1; ... % n2/n3
               1 1 0 0 1 1];    % n2/n3
             % 1 1 0 0 1 1 + low-alpha = rem
             % 1 0 0 0 1 1 + low-alpha = rem
           
% State transition matrix params [p_on p_off p_stay]
p_transition = [.2 .3 .8]; % [.1 .2 .9]; %

% Observation noise covariance
diag_R = 0.5*ones(length(sfreqs),1);
% Adjust for areas of higher noise
idx_sixtyHz = sfreqs>55 & sfreqs<65;
diag_R(idx_sixtyHz) = 3;

% State noise and initial state covariances
num_params = 0;
for ii = 1:length(peaks)
   num_params = num_params + peaks(ii).numPeakParams; 
end
diag_Q = 0.1*ones(num_params,1);
diag_CovX0 = 0.5*ones(num_params,1);
cnt = 0;
for ii = 1:length(peaks)
    for jj = 1:peaks(ii).numPeakParams
        switch lower(peaks(ii).modelParams(jj).linkFunClass)
            case 'exponential'
                cnt = cnt + 1;
                diag_Q(cnt) = 0.01;
                diag_CovX0(cnt) = 0.01;
            case 'sigmoid'
                cnt = cnt + 1;
                diag_Q(cnt) = 0.1;
                diag_CovX0(cnt) = 0.1;
            case 'identity'
                cnt = cnt + 1;
                diag_Q(cnt) = 1; 
                diag_CovX0(cnt) = 0.5;
        end
    end
end
CovX0 = diag(diag_CovX0);

%*************************
% Form state-space model *
%*************************
ss = StateSpaceMultiPeak(peaks,is_dynamic,peak_combos,[],[],[],[],[],diag_Q,sfreqs,diag_R);

%****************************************
% Add initial parameter values in model *
%****************************************
% Initial parameters for non-baseline peaks
%       Freq         Amp                        BW      SKEW
x0_hat = [0       ~f_link_fun_amp(1)*log(10)      0              ... %60Hz
          0       ~f_link_fun_amp(2)*log(.1)      -3       1      ... %Delta
          0       ~f_link_fun_amp(3)*log(.1)      0        0     ... %Alpha + Harmonics
          0       ~f_link_fun_amp(4)*log(.1)      -3            ... %Sigma
          0       ~f_link_fun_amp(5)*log(5)      0       1 ];     ... %Slow

% Initial paramteters for baseline
baseline_estimate = exp_baseline_fit(sfreqs,spect', 0);
EX0_bl = zeros(p_base.numPeakParams,1);
for ii = 1:p_base.numPeakParams
    EX0_bl(ii) = p_base.modelParams(ii).linkInvFun(baseline_estimate(ii));
end

% Set initial state in model
ss.EX0(1:(end-3)) = x0_hat;
ss.EX0((end-2):end) = EX0_bl;

% Set initial combo to wake in model
ss.alpha0 = [1 0 0 0 0 0]';

% Set initial state covariance in model
ss.CovX0 = CovX0;

%***********
% Run IEKF *
%***********
num_iters = 20;
num_particles = 2000;
verbose = 2;
ss_est = iekf(spect,stimes,is_artifact,ss,num_iters,num_particles,verbose);
% Add times to estimate structure 
ss_est.t = [stimes(1)-(1/Fs), stimes]; % 0:length(stimes);

%****************
% Plot and save *
%****************
if ~isempty(save_fname)
    save([save_fname '.mat']); % ,'sim_spect','stimes','sfreqs','X');
    plotMPSpectWithSlices(stimes,sfreqs_trunc,nanpow2db(spect_trunc),[],ss,ss_est,2,save_fname);
    plotMPStateEstimates(ss_est,ss,false,[],[],[],false,true,save_fname);
else
    plotMPSpectWithSlices(stimes,sfreqs_trunc,nanpow2db(spect_trunc),[],ss,ss_est,2);
    plotMPStateEstimates(ss_est,ss,false,[],[],[],false,true);
end
