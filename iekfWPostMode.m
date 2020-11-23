function [ss_estim] = iekfWPostMode(y, t, is_artifact, m, num_iters,  num_particles, verbose)
%iekfWPostMode estimates a multi-peak, state-space model for a spectrogram using
% the iterated extended Kalman filter augmented to estimate the discrete
% transition of the On/Off-peak combination and to draw an improved 
% initial reference trajectory. The difference from iekf is that it 
% uses an approximation to the posterior density to evaluate the draws and 
% select the intial reference trajectory instead of the likelihood.
% The posterior approximation is an internal function logPostx.
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
%   num_iters     -- positive integer number of iterations. 
%                    EKF version given by 1. Default 20.
%   num_particles -- positive integer number of draws. 
%                    Version without draws given by 1. Default 500.
%   verbose       -- flag indicating level of verbosity. Default 2.
%                    0 - no display of progress.
%                    1 - text display of progress
%                    2 - graphical display of fits as it progresses.
% OUTPUTS:
%   ss_etim -- estimate structure for StateSpaceMultiPeak. 
%              Contains:
%                num_particles
%                num_iters
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
if nargin < 7
    verbose = [];
end
if nargin < 6
    num_particles = [];
end
if nargin < 5
    num_iters = [];
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
    if isempty(num_particles)
        num_particles = 500;
    end
    if isempty(num_iters)
        num_iters = 20;
    end
    if isempty(t)
        t = 1:size(y,2);
    end
    T = length(t);
    if isempty(is_artifact)
        is_artifact = false(1,T);
    end
    
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
    
    
    if num_particles > 1
        f_cast_net = true;
        if verbose > 0
            disp(['Running with ' num2str(num_particles) ' particles']);
        end
    else
        f_cast_net = false;
    end
    
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
    
    % Initialize progressbar and figure
    if verbose == 2
        progressbar();
        fh = figure;
    end
    
    start_time = tic;
    
    %Number of time steps at the begining to do high iterations, set to
    %negative to remove
    warm_up = -1;
    warm_up_iterations = 100;
    
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
    xf_next = cell(m.numCombos,1);
    yf_next = cell(m.numCombos,1);
    comps_f_next = cell(m.numCombos,1);
    Udraws = cell(m.numCombos,1);
    
    % Extracts/forms these matrices for faster computation
    % ASSUMES observation noise covariance R is constant matrix
    omega = m.omega;
    Omega = repmat(omega,num_particles, 1)';
    R = m.R(1);
   
    % Use these two lines if using parfor
    Omega_parconst = parallel.pool.Constant(Omega);
    R_parconst = parallel.pool.Constant(R);
% If using for instead, comment two lines above and use those below
%     Omega_parconst.Value = (Omega);
%     R_parconst.Value = (R);

    %************
    % Time step *
    %************
    for kk = 2:(T+1)
        tic
        if verbose > 0
            iteration_time = tic;
            disp(['Time ' num2str(kk-1) ' out of ' num2str(T)]);
        end
        
        if kk < warm_up
            iterations = warm_up_iterations;
        else
            iterations = num_iters;
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
                
                xf_next{jj} = xp_curr;
                Pf_curr{jj} = Pp_curr;
            end
            yf_next = yp_curr;
            comps_f_next = comps_p_curr;
            alpha(:,kk) = alpha(:,kk-1);
            [~,mode_combo] = max(alpha(:,kk));
            
        else
            % Make draws for improved reference trajectory
            if f_cast_net
                parfor jj = 1:m.numCombos
                    Udraws{jj} = randn(num_particles-1,length(xp_curr));
                end
            end
            % Filter state update for each combo
            parfor jj = 1:m.numCombos
                curr_combo = m.peakCombos(jj,:);
                [yp_curr{jj}, comps_p_curr{jj}] = m.getH(xp_curr, omega, curr_combo);
                sse_pred(jj) = -logPostx(zeros(size(xp_curr)),Pp_curr,y_k-yp_curr{jj},R_parconst.Value);
                    
                % Use the predicted state as the reference point for the linearization
                eta_next = xp_curr;
                
                % Use sampling from predicted distribution to select alternative linearization point
                if f_cast_net
                    %Draw particles
                    x_draws = repmat(eta_next',num_particles,1);
                    w_draws = Udraws{jj}*sqrtP; % randn(num_particles-1,length(eta_next))*sqrtP;
                    w_draws(:,m.offIdxs{jj}) = 0;
                    x_draws(2:end,:) = x_draws(2:end,:) + w_draws;
                    w_draws = [zeros(1,m.dimX); w_draws];
                    
                    %Compute function for each particle
                    yhat_draws = m.getH(x_draws', Omega_parconst.Value, curr_combo);
                    resid_draws = repmat(y_k,[1 num_particles]) - yhat_draws;
                    sse_draws = -logPostx(w_draws',Pp_curr,resid_draws,R_parconst.Value);
                    
                    %Take maximum approx. posterior mode particle
                    [~,idx_min_sse] = min(sse_draws);
                    eta_next = x_draws(idx_min_sse,:)';
                    sse_draw = sse_draws(idx_min_sse);
                else
                    sse_draw = sse_pred(jj);
                end
                
                eta_curr = eta_next;
                sse_curr(jj) = sse_draw;
                
                % Perform iterations
                for ii = 1:iterations
                    eta_curr = eta_next;
                    
                    %Derivative of the observation
                    M_curr{jj} = m.getHx(eta_curr,omega,curr_combo);
                    if any(~isreal(M_curr{jj}))
                       disp('imag'); 
                    end
                    
                    %Computing the Kalman gain, K
                    MP_term = Pp_curr *M_curr{jj}';
                    block1 = M_curr{jj} * MP_term + R_parconst.Value;
                    K_curr = MP_term / block1; %** matrix close to singular
                    
                    %Computing the model estimate
                    [yf_curr{jj}, comps_f_curr{jj}] = m.getH(eta_curr, omega, curr_combo);
                    
                    %Update the reference point for the linearized filter estimate
                    eta_next  = xp_curr + K_curr * (y_k - yf_curr{jj} - M_curr{jj} * (xp_curr - eta_curr));
                end
                
                % Last eta_next is ultimately the filter estimate. But 
                % eta_curr, M_curr, yf_curr, and Pp_curr need to be retained 
                % to compute the data likelihood of the combo.
                xf_curr{jj} = eta_curr;
                xf_next{jj} = eta_next;
                [yf_next{jj}, comps_f_next{jj}] = m.getH(eta_next, omega, curr_combo);
                    
                block2 = eye(m.dimX) - K_curr * M_curr{jj};
                Pf_curr{jj} = block2 * Pp_curr * block2' + K_curr * R_parconst.Value * K_curr';
                sse_curr(jj) = -logPostx(xp_curr-eta_curr,Pp_curr,y_k-yf_curr{jj},R_parconst.Value);
            end
            
            %***********************
            % Determine mode combo *
            %***********************
            % Check if any filters improved on the predictions
            if min(sse_pred) < min(sse_curr)
                % If not, find "max posterior mode" prediction
                combo_ind=find(sse_pred==min(sse_pred));
                if length(combo_ind)>1
                    combo_ind=combo_ind(1);
                end
                if verbose > 0
                    disp(['***all filters have larger sse than prediction sse of combo ' num2str(combo_ind)]);
                end
                % No filtering. Set to best prediction.
                for jj = 1:m.numCombos
                    xf_next{jj} = xp_curr;
                    Pf_curr{jj} = Pp_curr;
                end
                yf_next = yp_curr;
                comps_f_next = comps_p_curr;
                alpha(combo_ind,kk) = 1;
                mode_combo = combo_ind;
                
            else
                % If so, check if we have a good reference for log likelihood ratio
                
                % Intialize reference combo residual and covariance
                ref_resid=[];
                ref_cov=[];
                good_ref=false;
                % Start with the last combo as the start and loop through to see if there is a good reference
                for ref_combo_check = m.numCombos:-1:1
                    ref_resid = y_k-yf_curr{ref_combo_check}+M_curr{ref_combo_check}*xf_curr{ref_combo_check}-M_curr{ref_combo_check}*xp_curr;
                    ref_cov = M_curr{ref_combo_check} * Pp_curr * M_curr{ref_combo_check}' + R;
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
                        xf_next{jj} = xp_curr;
                        Pf_curr{jj} = Pp_curr;
                    end
                    yf_next = yp_curr;
                    comps_f_next = comps_p_curr;
                    alpha(:,kk) = alpha(:,kk-1);
                    [~,mode_combo] = max(alpha(:,kk));
                    
                else
                    % If there IS a good reference compute the log likelihood ratio with respect to the reference
                    ref_fsse = .5*ref_resid'*(ref_cov\ref_resid);
                    for jj = setdiff(1:m.numCombos,ref_combo)
                        resid = y_k-yf_curr{jj}+M_curr{jj}*xf_curr{jj}-M_curr{jj}*xp_curr;
                        block1 = M_curr{jj} * Pp_curr *M_curr{jj}' + R;
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
            
        end
        
        %*************************************************
        % Collapse to filter estimates to mode combo and *
        % move current cell info to output matrices      *
        %*************************************************
        xp_hat(:,kk) = xp_curr;
        Pp(:,:,kk) = Pp_curr;
        yp_hat(:,kk) = yp_curr{mode_combo};
        comps_hat_p(:,:,kk) = comps_p_curr{mode_combo};
        xf_hat(:,kk) = xf_next{mode_combo}; 
        Pf(:,:,kk) = Pf_curr{mode_combo};
        yf_hat(:,kk) = yf_next{mode_combo}; 
        comps_hat_f(:,:,kk) = comps_f_next{mode_combo}; 
        
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
                uistack(ph,'top');
                
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
ss_estim.num_particles = num_particles;
ss_estim.num_iters = num_iters;
ss_estim.xf_hat = xf_hat;
ss_estim.Pf = Pf;
ss_estim.xp_hat = xp_hat;
ss_estim.Pp = Pp;
ss_estim.yp_hat = yp_hat;
ss_estim.yf_hat = yf_hat;
ss_estim.comps_hat_p = comps_hat_p;
ss_estim.comps_hat_f = comps_hat_f;
ss_estim.alpha = alpha;

end

%***************************************************************************
% Function to evaluate the log of the approximate filter posterior density *
%***************************************************************************
function p = logPostx(dx,P_tm1,resid,R)
% logPostx evaluates the log of an approximation to the filter posterior 
% density. It is used to compare the predictive draws and select the modal
% draw as the initial reference trajectory.
%
% The approximation is N(y-h(x'),R)*N(x'-x(t|t-1),P(t|t-1)).
%
% INPUTS:
%   dx    -- matrix (dim_x x num_draws) of state inputs, i.e., deviations
%            of draws x' from x(t|t-1), not the draws themselves.  
%   P_tm1 -- P(t|t-1) prediction error covariance matrix (dim_x x dim_x). 
%   resid -- matrix (dim_y x num_draws) of observation residuals, i.e.,
%            y(t)-h(x(t|t-1)+dx). 
%   R     -- observation noise covariance matrix (dim_y x dim_y).
%   All inputs required.
% Outputs:
%   p -- a vector (num_draws x 1) containing density values for each draw
%

    dim_x = size(dx,1);
    dim_y = size(resid,1);
    logpred_const = (-dim_x/2)*log(2*pi) + (-1/2)*log(det(P_tm1)); 
    loglik_const = (-dim_y/2)*log(2*pi) + (-1/2)*log(det(R));
    p = -1/2*sum((dx'/P_tm1).*dx',2) + -1/2*sum((resid'/R).*resid',2);
    p = p + loglik_const + logpred_const;
end
