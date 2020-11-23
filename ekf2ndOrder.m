function [ss_estim] = ekf2ndOrder(y, t, is_artifact, m, verbose, ss_estim, t_start,f_trunc,f_modified)
%ekf2ndOrder estimates a multi-peak, state-space model for a spectrogram using
% a 2nd-order extended Kalman filter, augmented to estimate the discrete
% transition of the On/Off-peak combination. Two optional forms are 
% available---the modified truncated 2nd-order filter and the modified 
% 2nd-order Gaussian filter---based on Appendix 9B from Jazwinski. 
%
% The non-modified forms were not completed due to the tediousness of the 
% computation of the D array. They can be completed by including the 
% computation of D below. 
% 
% The modified Gaussian filter appears quite stable and performs generally
% well. But we have found the modified truncated filter to be fairly 
% unstable, resulting in eventual crashes. 
% 
% NOTE: This function explicitly overwrites the state transition by 
%       incorporating the 0.9 factor into Phi and using this Phi for 
%       both the state and state-error covariance updates. 
%       
% INPUTS: 
%   y             -- observed spectrogram (dim_y x T)
%   t             -- time points of spectrogram (1 x T)
%   is_artifact   -- indicator vector of time points of artifacts (1 x T)
%   m             -- StateSpaceMultiPeak object containing the model
%   verbose       -- flag indicating level of verbosity. Default 2.
%                    0 - no display of progress.
%                    1 - text display of progress
%                    2 - graphical display of fits as it progresses.
%   ss_estim      -- state estimate structure from prior estimate
%   t_start       -- positive integer time step at which to start filter. 
%                    Default is 1.
%                 [ss_estim and t_start are legacy inputs from attempts to 
%                  debug the truncated filter.] 
%   f_trunc       -- flag indicating type of filter to use. Default is 0.
%                    0 - truncated. 1 - Gaussian. 
%   f_modified    -- flag indicating form of filter used. Default is 1.
%                    0 - non-modified form (NOT fully implemented)
%                    1 - modified form.
%
% OUTPUTS:
%   ss_etim -- estimate structure for StateSpaceMultiPeak. 
%              Contains:
%                xf_hat      - filter state estimates (dim_x x T+1)
%                Pf          - filter state error covariances (dim_x x dim_x x T+1)
%                yf_hat      - filter observation estimates (dim_y x T+1)
%                comps_hat_f - filter component peak estimates (dim_y x T+1)
%                alpha       - filter estimate of On/Off-peak combination (num_combos x T+1) 
%                xp_hat      - prediction state estimates (dim_x x T+1)
%                Pp          - prediction state error covariances (dim_x x dim_x x T+1)
%                yp_hat      - prediction observation estimates (dim_y x T+1)
%                comps_hat_p - prediction component peak estimates (dim_y x T+1)
%
%
% Created by Patrick Stokes and Michael Prerau
% Created on 2017-04-20
% Modified on 2017-04-24
% 

%*************************
% Handle variable inputs *
%*************************
if nargin < 9
    f_modified = [];
end
if nargin < 8
    f_trunc = [];
end
if nargin < 7
    t_start = [];
end
if nargin < 6
    ss_estim = [];
end
if nargin < 5
    verbose = [];
end
if nargin < 4
    m = [];
end
if nargin < 3
    is_artifact = [];
end
if nargin < 2 
    t = [];
end
if nargin < 1
    y = [];
end

%************************
% Check required inputs *
%************************
if isempty(m)
    error('multi-peak state-space model m not provided as input.');
elseif isempty(y)
    error('spectrogram observation data y not provided as input.');
elseif m.dimY ~= size(y,1) && m.dimY ~= size(y,2)
    error('dimension of observation in data and model are not equal.');
else
%*************
% Run filter *
%*************    
    % Set defaults
    if m.dimY ~= size(y,1)
       y = y'; 
    end
    if isempty(verbose)
        verbose = 2; % 0 - no display, 1 - text only, 2 - text and graphics
    end
    if isempty(t)
        t = 1:size(y,2);
    end
    T = length(t);
    if isempty(is_artifact)
        is_artifact = false(1,T);
    end
    
    if ~isempty(ss_estim)
        % Option to re-use previous estimates
        xf_hat = ss_estim.xf_hat;
        Pf = ss_estim.Pf;
        xp_hat = ss_estim.xp_hat;
        Pp = ss_estim.Pp;
        yp_hat = ss_estim.yp_hat;
        yf_hat = ss_estim.yf_hat;
        comps_hat_p = ss_estim.comps_hat_p;
        comps_hat_f = ss_estim.comps_hat_f;
        alpha = ss_estim.alpha;
        
    else
        % Allocate space for filter variables
        xp_hat = zeros(m.dimX, T+1);
        Pp = zeros(m.dimX, m.dimX, T+1);
        xf_hat = zeros(m.dimX, T+1);
        Pf = zeros(m.dimX, m.dimX, T+1);
        yp_hat = zeros(m.dimY, T+1);
        yf_hat = zeros(m.dimY, T+1);
        comps_hat_p = zeros(m.numPeaks, m.dimY, T+1);
        comps_hat_f = zeros(m.numPeaks, m.dimY, T+1);
        log_lik_combo = -inf(m.numCombos, T+1);
        log_lik_full = -inf(m.numCombos, T+1);
        alpha = zeros(m.numCombos, T+1);

        % Pad t and y to lineup with filter
        t = [0 t];
        y = [zeros(m.dimY,1) y];
        % Initialize the filter
        xp_hat(:,1) = m.EX0;
        Pp(:,:,1) = m.CovX0;
        xf_hat(:,1) = m.EX0;
        Pf(:,:,1) = m.CovX0;
        % Initialize the peak-on matrix
        alpha(:,1) = m.alpha0;
    end
    
    if isempty(t_start)
        t_start = 1;
    end
    
    if isempty(f_trunc)
        f_trunc = 0;
    end
    if isempty(f_modified)
        f_modified = 1;
    end
    
    % Initialize progressbar and figure
    if verbose == 2
        progressbar();
        fh = figure;
    end
    
    start_time = tic;
    
    % Used in graphical display of progress bar, negative means plot right away
    warm_up = -1;
    
    % Temporary storage for estimates at a given time to allow parfor
    % slicing
    sse_pred = zeros(m.numCombos,1);
    sse_curr = zeros(m.numCombos,1);
    yp_curr = cell(m.numCombos,1);
    comps_p_curr = cell(m.numCombos,1);
    xf_curr = cell(m.numCombos,1);
    Pf_curr = cell(m.numCombos,1);
    yf_curr = cell(m.numCombos,1);
    comps_f_curr = cell(m.numCombos,1);
    M_curr = cell(m.numCombos,1);
    curr_resid = cell(m.numCombos,1);
    Y = cell(m.numCombos,1);
    
    Hxx_curr = cell(m.numCombos,1);
    ExhTmExEhT_curr = cell(m.numCombos,1);
    Pd2h_curr = cell(m.numCombos,1);
    d2hP2d2h_curr = cell(m.numCombos,1);
    E_hmEhhmEhT_curr = cell(m.numCombos,1);
    
    % Extracts/forms these matrices for faster computation
    % ASSUMES observation noise covariance R is constant matrix
    omega = m.omega;
    R = m.R(1);
% Use this line if using parfor
    R_parconst = parallel.pool.Constant(R);
% If using for instead, comment line above and use one below
%     R_parconst.Value = (R);
  
    %************
    % Time step *
    %************
    for kk = (t_start+1):(T+1)
        tic
        if verbose > 0
            iteration_time = tic;
            disp(['Time ' num2str(kk-1) ' out of ' num2str(T)]);
        end
        
        % Temporary storage for convenience and to improve parfor
        y_k = y(:,kk);
        Q_k = m.Q(t(kk));
        Phi_k = m.Phi(t(kk),t(kk-1));
        % fintegr_k = m.fIntegr(xf_hat(:,kk-1),t(kk-1),t(kk),m.fIntegrParams{:});
        xf_old = xf_hat(:,kk-1);
        Pf_old = Pf(:,:,kk-1);
        
        %******************
        % Prediction step *
        % *****************
        % Same for all combos.
        % NOTE: That state transition is being over-ruled
        %       0.9 factor added to Phi and state uses Phi instead of fintegr 
        Phi_k = 0.9 * Phi_k;
        xp_curr = Phi_k * xf_old; % xp_curr = xf_old + fintegr_k;
        Pp_curr = Phi_k * Pf_old * Phi_k' + Q_k;
        sqrtP = chol(Pp_curr);
        
        %*************************************
        % Filter state and combo update step *
        %*************************************
        if is_artifact(kk-1) % is_artifact was not augmented with zero in front, so is lagging.
            if verbose>0
                disp('Artifact!');
            end
            % No filter update. Use max likelihood prediction combo.
            parfor jj = 1:m.numCombos
                curr_combo = m.peakCombos(jj,:);
                [yp_curr{jj}, comps_p_curr{jj}] = m.getH(xp_curr, omega, curr_combo);
                sse_pred(jj) = (y_k-yp_curr{jj})'*(R\(y_k-yp_curr{jj}));
                
                xf_curr{jj} = xp_curr;
                Pf_curr{jj} = Pp_curr;
            end
            yf_curr = yp_curr;
            comps_f_curr = comps_p_curr;
            alpha(:,kk) = alpha(:,kk-1);
            [~,mode_combo] = max(alpha(:,kk));
            
        else
            % Filter state update for each combo
            parfor jj = 1:m.numCombos
                curr_combo = m.peakCombos(jj,:);
                
                % Current observation function
                [yp_curr{jj}, comps_p_curr{jj}] = m.getH(xp_curr, omega, curr_combo);
                sse_pred(jj) = (y_k-yp_curr{jj})'*(R_parconst.Value\(y_k-yp_curr{jj}));
                
                % Current 1st derivative of the observation function
                M_curr{jj} = m.getHx(xp_curr,omega,curr_combo);
                if any(~isreal(M_curr{jj}))
                    disp('imag');
                end
                
                % Current 2nd derivative of the observation function
                Hxx_curr{jj} = m.getHxx(xp_curr,omega,curr_combo);
                
                % Matrices and terms needed for updates
                ExhTmExEhT_curr{jj} = Pp_curr*M_curr{jj}';                
                Pd2h_curr{jj} =  squeeze(sum(sum(repmat(Pp_curr,[1 1 length(omega) ]).* permute(Hxx_curr{jj},[2 3 1]),2),1));
                if f_trunc
                    % Truncated
                    E_hmEhhmEhT_curr{jj} = M_curr{jj}*Pp_curr*M_curr{jj}' - (1/4)*(Pd2h_curr{jj}*Pd2h_curr{jj}');
                else
                    % Gaussian
                    d2hP2d2h_curr{jj} = zeros(length(omega),length(omega));
                    for mm = 1:length(omega)
                        term = (Pp_curr*squeeze(Hxx_curr{jj}(mm,:,:))*Pp_curr')';
                        d2hP2d2h_curr{jj}(:,mm) = sum(sum(Hxx_curr{jj}.*permute(repmat(term,[1 1 length(omega)]),[3 1 2]),3),2);
                    end
                    E_hmEhhmEhT_curr{jj} = M_curr{jj}*Pp_curr*M_curr{jj}' + (1/2)*d2hP2d2h_curr{jj};
                end
                Y{jj} = E_hmEhhmEhT_curr{jj}+R_parconst.Value;
                
                a = xp_curr;
                B = Pp_curr*M_curr{jj}'/Y{jj};
                
                % Residual used in log-likelihood ratio of combos
                curr_resid{jj} = y_k - yp_curr{jj} - (1/2)*Pd2h_curr{jj};
                
                % Compute filter state and state-error covariance estimates
                xf_curr{jj} = a + B * curr_resid{jj};
                C = Pp_curr - B*M_curr{jj}*Pp_curr;
                Pf_curr{jj} = C;
                % Non-modified forms were not finished due to tediousness
                % of computation of D. To finish D must be computed.
                if ~f_modified
                    D = 0; %***** Have not implemented computation of D (see Jazwinski p.364)
                    for ll = 1:length(omega)
                        Pf_curr{jj} = Pf_curr{jj} + D(:,:,ll)*curr_resid{jj}(ll);
                    end
                end
                [yf_curr{jj}, comps_f_curr{jj}] = m.getH(xf_curr{jj}, omega, curr_combo);
                sse_curr(jj) = (y_k-yf_curr{jj})'*(R_parconst.Value\(y_k-yf_curr{jj}));
            end
            %         tocBytes(gcp)
            %
            % Log likelihood ratio with better handling of linearization
            % Check if any filters improved on the prediction
            % Check if we found a good reference for log likelihood ratio
            % Intialize reference combo residual and covariance
%             if min(sse_pred) < min(sse_curr)
%                 combo_ind=find(sse_pred==min(sse_pred));
%                 if length(combo_ind)>1
%                     combo_ind=combo_ind(1);
%                 end
%                 if verbose > 0
%                     disp(['***all filters have larger sse than prediction sse of combo ' num2str(combo_ind)]);
%                 end
%                 % No filtering. Set to best prediction.
%                 for jj = 1:m.numCombos
%                     xf_curr{jj} = xp_curr;
%                     Pf_curr{jj} = Pp_curr;
%                 end
%                 yf_curr = yp_curr;
%                 comps_f_curr = comps_p_curr;
%                 alpha(combo_ind,kk) = 1;
%                 mode_combo = combo_ind;
%                 
%             else

            %************************************************
            % Determine mode combo via log likelihood ratio *
            %************************************************
            % Intialize reference combo residual and covariance
                ref_resid=[];
                ref_cov=[];
                good_ref=false;
                % Start with the last combo and loop through to see if there is a good reference
                for ref_combo_check = m.numCombos:-1:1
                    ref_resid = curr_resid{ref_combo_check};
                    ref_cov = Y{ref_combo_check};
                    if det(ref_cov)>0
                        log_lik_combo(ref_combo_check,kk)=0;
                        ref_combo=ref_combo_check;
                        good_ref=true;
                        if verbose > 0
                            disp(['Setting combo ' num2str(ref_combo) ' as reference combo']);
                        end
                        break;
                    end
                end
                
                % If there is no good reference, treat it like an artifact and move
                % on to the next time point
                if ~good_ref
                    if verbose > 0
                        disp('All combos bad, setting to artifact!!!!');
                    end
                    % No filtering
                    for jj = 1:m.numCombos
                        xf_curr{jj} = xp_curr;
                        Pf_curr{jj} = Pp_curr;
                    end
                    yf_curr = yp_curr;
                    comps_f_curr = comps_p_curr;
                    alpha(:,kk) = alpha(:,kk-1);
                    [~,mode_combo] = max(alpha(:,kk));
                    
                else
                    % If there IS a good reference compute the log likelihood ratio with respect to the reference
                    ref_fsse = .5*ref_resid'*(ref_cov\ref_resid);
                    for jj = setdiff(1:m.numCombos,ref_combo)
                        resid = curr_resid{jj};
                        block1 = Y{jj};
                        log_lik_combo(jj,kk)=-.5*resid'*(block1\resid)...
                            +ref_fsse+1/2*log(det(ref_cov)/det(block1));
                    end
                    
                    % Incorporate transition probabilities into log likelihood ratio
                    log_lik_full(:,kk) = log_lik_combo(:,kk) + log(m.transMatr*alpha(:,kk-1));
                    log_lik_full(:,kk) = log_lik_full(:,kk) - log(m.transMatr(ref_combo,:)*alpha(:,kk-1));
                    
                    % Determine posterior mode combo from log likelihood ratio
                    [~,mode_combo] = max(log_lik_full(:,kk));
                    alpha(mode_combo,kk) = 1;
                end
            end
            
%        end
        
        %*************************************************
        % Collapse to filter estimates to mode combo and *
        % move current cell info to output matrices      *
        %*************************************************
        xp_hat(:,kk) = xp_curr;
        Pp(:,:,kk) = Pp_curr;
        yp_hat(:,kk) = yp_curr{mode_combo};
        comps_hat_p(:,:,kk) = comps_p_curr{mode_combo};
        xf_hat(:,kk) = xf_curr{mode_combo};
        Pf(:,:,kk) = Pf_curr{mode_combo};
        yf_hat(:,kk) = yf_curr{mode_combo};
        comps_hat_f(:,:,kk) = comps_f_curr{mode_combo};
        
        %******************
        % Output displays *
        %******************
        if verbose>0
            % The current most likely peak combo
            disp(['current mixture: ' num2str(alpha(:,kk)')]);
            
            % Finish time estimate
            if ~mod(kk-1,round(T/300)) || kk==2
                disp([num2str(100*(kk-1)/T) '%  ' num2str((toc(start_time)/(kk-1))*(T-(kk-1))/60) ' minutes left']);
            end
            
            % Graphical displays
            if verbose == 2
                % Progress bar
                if kk > warm_up
                    progressbar((kk-1)/T);
                end
                
                % Plot
                figure(fh);
                cla
                hold all;
                plot(omega, squeeze(comps_hat_f(:,:,kk))','linewidth',3);
                th = plot(omega,y(:,kk),'k','linewidth',4);
                ph = plot(omega, squeeze(yf_hat(:,kk,1))','linewidth',3);
                set(ph,'linewidth',6,'color','b');
                
                uistack(th,'bottom');
                uistack(ph,'bottom');
                
                title(['Time = ' num2str(kk-1) ' Peaks On: [' num2str(m.peakCombos(mode_combo,:)) ']']);
                drawnow;
            end
            
            % Iteration time
            disp(['Iteration took ' num2str(toc(iteration_time)) ' seconds']);
        end
    end
    
end

%******************************************
% Store estimates in structure for output *
%******************************************
ss_estim.xf_hat = xf_hat;
ss_estim.Pf = Pf;
ss_estim.xp_hat = xp_hat;
ss_estim.Pp = Pp;
ss_estim.yp_hat = yp_hat;
ss_estim.yf_hat = yf_hat;
ss_estim.comps_hat_p = comps_hat_p;
ss_estim.comps_hat_f = comps_hat_f;
ss_estim.alpha = alpha;
