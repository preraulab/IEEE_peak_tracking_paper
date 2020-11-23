function compareEKFs_jobs_wRandSSMs(odir,sim_type,num_peaks,num_sims,NT,dT,maxF,NF,rand_seeds,num_particles,num_iters,K_lim,T_lim,f_verb)
% compareEKFs_jobs_wRandSSMs compares the performance of the variant IEKF
% filters by generating simulations, estimating all filters on each simulation, 
% and computing a set of 16 statistical measures and tests from the estimates.  
%
% NOTE: The first four, primary inputs are as strings to facilitate batch 
%       submission to a cluster (with other inputs as default). 
%
% NOTE: This implementation assumes 4 cores available for parallel
%       processing. This is set in the variable num_cores. It can be increased 
%       or decreased as necessary.
% 
% Inputs:
%   odir          -- string of output directory path where output is saved.
%                    Default is 'results'.
%   sim_type      -- string indicating 'randWalk' (default) or 'pseudoDeterm'.
%   num_peaks     -- string of integer number of peaks on an exponential-decay background. 
%                    Default is '1'.
%   num_sims      -- string of number of simulations. Default is '1'.
%   NT            -- integer number of time steps for length of the
%                    simulations. Default is 100. 
%   dT            -- sampling rate of simultion. It is completely arbitrary
%                    and does not affect anything beyond the label of the
%                    time bins. Default is 0.005
%   maxF          -- positive Nyquist (maximum) frequency of simulations.
%                    Default is 100.
%   NF            -- integer number of equally spaced frequency bins 
%                    from 0 to maxF. Default is 300.
%   rand_seeds    -- an array of 1xnum_simsx10 of random seed values.
%                    The 10 random seed values for each simulation 
%                    are initializations for 
%                    1. the random setting of peak types and frequency ranges,
%                    2. the generation of the simulation, and
%                    3-10. the separate filter estimates.
%                    *WARNING* The re-use of seeds does not appear to work.
%                    At some point within the code the seeds are being
%                    overridden.
%   num_particles -- an integer number of particle draws to obtain an improved
%                    initial reference trajector. Used in the iekf, 
%                    iekfW2ndDeriv, and the iekfWPostMode. Basic EKF and 
%                    IEKF versions obtained by using num_particles=1. 
%                    Default is 1000.
%   num_iters     -- an integer number of iterations used in the iekf, 
%                    iekfW2ndDeriv, and the iekfWPostMode. Corresponding 
%                    EKFs obtained by using num_iters=1. Default is 10.
%   K_lim         -- an integer limit on the number on the cross-frequency
%                    correlations evaluated in the checkEstimates
%                    computation of iCorr, boxQ, and moranI. Default is 4.
%   T_lim         -- an integer limit on the number on the time-lag
%                    correlations evaluated in the checkEstimates
%                    computation of boxQ and moranI. Default is 3.
%   f_verb        -- integer flag for verbosity of each filter computation. 
%                    0 -- no text or graphical display of filter progress.
%                    1 -- text display of filter progress.
%                    2 -- text and graphical display of filter.
%                    Default is 0.
%
% Outputs: Outputs are saved to a file in odir. 
%   Output file name is set as [odir '/peaks' num2str(num_peaks) '_' curr_date_str '.mat']; 
%   stats      -- 8x16xnum_sims array of stats for each filter
%                 estimate of each simulation. There are 16 statistical
%                 measures computed for each of the 8 filters on each
%                 simulation.
%   rand_seeds -- array of generated random seeds. *WARNING* as mentioned
%                 above the re-use of the random seeds does not appear
%                 functional at this time.
%   Some of the input settings are also saved. 
%
% Created -- Patrick Stokes
%

%*************************
% Handle variable inputs *
%*************************
if nargin < 14
    f_verb = [];
end
if nargin < 13
    T_lim = [];
end
if nargin < 12
    K_lim = [];
end
if nargin < 11
    num_iters = [];
end
if nargin < 10
    num_particles = [];
end
if nargin < 9
    rand_seeds = [];
end
if nargin < 8
    NF = [];
end
if nargin < 7
    maxF = [];
end
if nargin < 6
    dT = [];
end
if nargin < 5
    NT = [];
end
if nargin < 4
    num_sims = [];
end
if nargin < 3
    num_peaks = [];
end
if nargin < 2
    sim_type = [];
end
if nargin < 1
    odir = []; 
end
    
%*********************
% Set default inputs *
%*********************
if isempty(odir)
    odir = 'results';
end
if isempty(sim_type)
   sim_type = 'randWalk'; 
end
if isempty(num_peaks)
   num_peaks = '1'; 
end
if isempty(num_sims)
   num_sims = '1';
end
if isempty(NT)
    NT = 100;
end
if isempty(dT)
    dT = 0.005;
end
if isempty(maxF)
    maxF = 100;
end
if isempty(NF)
    NF = 300;
end
if isempty(num_particles)    
    num_particles = 1000;
end
if isempty(num_iters)
    num_iters = 10;
end
if isempty(K_lim)
    K_lim = 4;
end
if isempty(T_lim)
    T_lim = 3;
end
if isempty(f_verb)
    f_verb = 0;
end

%**********************
% Check/start parpool *
%**********************
num_cores = 4;
% Check if parpool exists. If not, do not create new one.
p = gcp('nocreate'); 
% If current number of workers is less than the number of cores, delete and recreate. 
if ~isempty(p)
    if p.NumWorkers < num_cores
        p.delete;
        p = gcp('nocreate');
    end
end
% If no parpool, start one.
if isempty(p)
    % Wait for parpool license to become available
    while ~license('checkout','Distrib_Computing_Toolbox') % license('checkout','parpool')
        pause(5);
    end
    parpool('local',num_cores);
end

%*************************************************
% Extract and setup needed variables and storage *
%*************************************************
num_sims = str2double(num_sims);
num_peaks = str2double(num_peaks);
curr_date_str = myDateStr(5);
ofile_name = [odir '/peaks' num2str(num_peaks) '_' curr_date_str '.mat'];
num_ests = 8;
num_stats = 16;
if isempty(rand_seeds)
    rand_seeds = randi(10^7,[length(num_peaks),num_sims,num_ests+2]);
end
stats = zeros(num_ests,num_stats,num_sims,length(num_peaks));
pick_combos_array = [];

% Generate simulations, estimating each filter variant and computing its 
% statistical measure on each simulatation.
% NOTE: Because num_peaks is a string input of a single integer and 
%       converted to a number, length(num_peaks)=1 and the for loop 
%       is not necessary. The for loop is present because the original code 
%       version used a vector of integers, e.g. [2 5 3], for input to num_peaks.
for ii = 1:length(num_peaks)
    disp(['Running with ' num2str(num_peaks(ii)) ' peaks...']);
    
    % All peaks except the exponential-decay background can dynamically turn On/off
    is_dynamic = [0 ones(1,num_peaks(ii))];
    
    % Temporary storage to save for true On/Off-comination and its estimates.
    tmp_pick_combos = zeros(2^(sum(is_dynamic)),NT,num_ests+1,num_sims);
    
    for jj = 1:num_sims
        disp(['  Simulation number ' num2str(jj) ' of ' num2str(num_sims) '...']);
        
        % Randomly set peak types and frequency ranges
        rng(rand_seeds(ii,jj,1));
        peak_types = rand(num_peaks(ii),1);
        freq_lens = rand(num_peaks(ii),1) + 1/num_peaks(ii);
        freq_lens = freq_lens/sum(freq_lens)*maxF;
        freq_ul = cumsum(freq_lens);
        freq_ll = [0; freq_ul(1:(end-1))];
        freq_lims = [freq_ll freq_ul];
        
        % Set up peak models and parameter bounds based on peak type
        peaks = PeakModel('expdecay',{'sigmoid','sigmoid','sigmoid'},{[20 40],[0.01 .2],[-20 0]},[],'Baseline');
        for kk = 1:num_peaks(ii)
            bw_lims = [3 min(49,max(5,(freq_lens(kk)/2)^2))];
            if peak_types(kk) < .25
                new_peak = PeakModel('skewgamma',[],{freq_lims(kk,:),[4*sqrt(1) 20],bw_lims,[.5 1]},[],[num2str(kk) ' Skew']);
            elseif peak_types(kk) < .75
                num_harm = floor(maxF/freq_lims(kk,2))-1;
                if num_harm > 0
                    new_peak = PeakModel('harmonic',[],{freq_lims(kk,:),[4*sqrt(1) 20],bw_lims,[.2 .7]}, num_harm,[num2str(kk) ' Harm']);
                else
                    new_peak = PeakModel('gaussian',[],{freq_lims(kk,:),[4*sqrt(1) 20],bw_lims}, [],[num2str(kk) ' Gauss']);
                end
            else
                new_peak = PeakModel('box',[],{freq_lims(kk,:),[4*sqrt(1) 20],bw_lims},4,[num2str(kk) ' Box']);
                
            end
            peaks = [peaks new_peak];
        end
        
        % Form state-space model and simulate
        stimes = dT*(1:NT);
        sfreqs = linspace(0,maxF,NF);
        rng(rand_seeds(ii,jj,2));
        if strcmpi('pseudoDeterm',sim_type)
            s = StateSpaceMultiPeak(peaks,is_dynamic,[],[],[],[],[],[],[],sfreqs);
            [sim_spect,X,alpha, is_on, V, comps] = simulatePseudoDeterm(s,[0 stimes]);
            sim_spect = sim_spect(:,2:end);
        else % 'randWalk'
            s = StateSpaceMultiPeak(peaks,is_dynamic,[],[],[.2 .2 .7],[],[],[],[],sfreqs);
            [sim_spect,X,alpha,is_on,t,W,V,comps] = s.simulate(NT,1);
        end
        % Temporarily store true On/Off-peak combinations
        tmp_pick_combos(:,:,1,jj) = alpha(:,2:end);
        
        % Estimate each filter variant
        for kk = 1:num_ests
            disp(['    Estimating filter ' num2str(kk) ' of ' num2str(num_ests) '...']);
            if kk == 1 
                % EKF with no draws
                tstart = tic;
                est = iekf(sim_spect,[],[],s,1,1,f_verb);
                rt = toc(tstart);
            elseif kk == 2 
                % EKF with draws
                tstart = tic;
                rng(rand_seeds(ii,jj,kk+2));
                est = iekf(sim_spect,[],[],s,1,num_particles,f_verb);
                rt = toc(tstart);
            elseif kk == 3
                % IEKF with no draws
                tstart = tic;
                rng(rand_seeds(ii,jj,kk+2));
                est = iekf(sim_spect,[],[],s,num_iters,1,f_verb);
                rt = toc(tstart);
            elseif kk == 4
                % IEKF with draws
                tstart = tic;
                rng(rand_seeds(ii,jj,kk+2));
                est = iekf(sim_spect,[],[],s,num_iters,num_particles,f_verb);
                rt = toc(tstart);
            elseif kk == 5
                % IEKF with 2nd-derivative, no draws
                tstart = tic;
                est = iekfW2ndDeriv(sim_spect,[],[],s,num_iters,1,f_verb);
                rt = toc(tstart);
            elseif kk == 6
                % IEKF with 2nd-derivative and draws
                tstart = tic;
                rng(rand_seeds(ii,jj,kk+2));
                est = iekfW2ndDeriv(sim_spect,[],[],s,num_iters,num_particles,f_verb);
                rt = toc(tstart);
            elseif kk == 7
                % 2nd-Order Gaussian EKF
                tstart = tic;
                est = ekf2ndOrder(sim_spect,[],[],s,f_verb);
                rt = toc(tstart);
            elseif kk == 8
                % EKF with draws using posterior-mode for selection 
                tstart = tic;
                est = iekfWPostMode(sim_spect,[],[],s,1,num_particles,f_verb);
                rt = toc(tstart);
            else
                
            end
            % Add time vector to estimate structure
            est.t = [(stimes(1)-dT) stimes];
            
            % Compute statistical measures for filter estimate
            stats(kk,1,jj,ii) = rt;
            stats(kk,2:end,jj,ii) = checkEstimates(est,s,sim_spect,X,alpha,s,K_lim,T_lim);
            %*************************************************************
            % NOTE: The first input of the model s into the above 
            %       is not technically correct because the 0.9 factor in 
            %       the state transition is hard coded into the filter 
            %       functions. So to be "correct" a new model s_est with 
            %       the 0.9 factor incorporated into the state transition
            %       should be formed and input. BUT only some constant 
            %       values are used, so it does not affect the result.
            %*************************************************************
            
            % Temporarily store estimated On/Off-peak combinations
            tmp_pick_combos(:,:,kk+1,jj) = est.alpha(:,2:end);
        end
    end
    % Store true and estimated On/Off-peak combinations
    if isempty(pick_combos_array)
        pick_combos_array.a = tmp_pick_combos;
    else
        pick_combos_array(ii,1).a = tmp_pick_combos;
    end
end

%**************
% Save output *
%**************
if ~isempty(ofile_name)
    % Make output directory if necessary
    fname_split = strsplit(ofile_name,'/');
    if length(fname_split)==1
        s = 1;
    else
        fpath = strjoin(fname_split(1:(end-1)),'/');
        if ~exist(fpath,'dir')
            disp([fpath ' does not exist. Making  ...']);
            [s,~,~] = mkdir(fpath);
            if s
                disp('  mkdir successful.');
            else
                disp('  mkdir unsuccessful. Output not saved.');
            end
        else
            s = 1;
        end
    end
    % Save
    if s > 0
        save(ofile_name,'stats','rand_seeds','pick_combos_array','num_peaks','num_sims','NT','dT','maxF','NF','num_particles','num_iters');
    end
end