classdef StateSpaceMultiPeak < MultiPeakModel
    %StateSpaceMultiPeak class extends the MultiPeakModel class to 
    % to form a state space model with multiple peak observation function
    % and enable its estimation. Contains a simulate function, 
    % 
    % NOTE: The class was built with the intention of allowing general 
    %       SDEs for the state evolution equation. However, this was never 
    %       fully realized, though the inputs are still in their functional 
    %       form (in particular, Phi, f_integr, and Q). Consequently, the 
    %       simulate class function is also not generalized. It assumes an 
    %       identity transition, constant noise covariances, and times 
    %       given by the integer times steps. Moreover, because of this and 
    %       for computation speed, the EKF/IEKF filters used to estimate the 
    %       model overwrite the state evolution equation, assuming  
    %       x(t) = 0.9*Phi*x(t-1) + w(t), with Phi(t1,t2) = I, and 
    %       constant noise covariances Q and R. 
    %
    %   obj = StateSpaceMultiPeak(peakModels, isDynamic, peakCombos, 
    %                 probTrans, f, Phi, fIntegr, fIntegrParams, Q, omega, 
    %                 R, EX0, CovX0, alpha0)
    %
    %   INPUTS:
    %     peakModels    -- an array (num_peaks x 1) of PeakModel objects.
    %     isDynamic     -- an indicator vector (1 x num_peaks) of which 
    %                      peaks are allowed to turn on and off. 
    %                      Default is ones(num_peaks,1).
    %     peakCombos    -- a indicator matrix (num_combos x num_peaks) of
    %                      which peaks are on in which combos. Default is
    %                      all possible combos from makePeakCombos.
    %     probTrans     -- scalar or vector of probability inputs to 
    %                      makeComboTransitionMatr. Can be 1x1, 1x2, or 1x3.
    %                      Default is [0.2 0.2 0.9].
    %     f             -- drift function of state f(t,x). Currently unused. 
    %                      Default is [].
    %     Phi           -- state transition function Phi(t2,t1) that takes 
    %                      an initial time t1 and finial time t2 and
    %                      returns a matrix (dim_x x dim_x).
    %                      Default is @(t2,t1)(eye(dim_x)).            
    %     fIintegr      -- integral of state drift function, integr_f(x,t1,t2,...)
    %                      takes state x, initial time t1, final time t2, and any additional parameters
    %                      returns matrix (dim_x x 1).
    %                      Default is @(x,t1,t2,params)(zeros(dim_x,1)).
    %     fIntegrParams -- additional parameters needed for evalulation of integr_f
    %                      takes final time t2 and initial time t1.
    %                      Default is empty {}.
    %     Q             -- state noise covariance function Q(t) that 
    %                      takes time t and returns matrix (dim_x x dim_x).
    %                      Default is @(t)(eye(dim_y)).
    %                      If Q is a scalar, a matrix with that value along
    %                      the diagonal is returned at each time.
    %                      If Q is a vector of same dimension as the 
    %                      observations, a matrix with that vector on the 
    %                      diagonal is returned at each time.
    %                      If Q is a matrix, that matrix is returned at
    %                      each time. 
    %                      Otherwise, Q is the provided function. 
    %     omega         -- a vector (1 x dim_y) of frequency bins at which
    %                      the observations are made and the observation 
    %                      function evaluated
    %     R             -- observation noise covariance function R(t) that
    %                      takes time t and and returns a matrix (dim_y x dim_y).
    %                      Default is @(t)(eye(dim_y)).
    %                      If R is a scalar, a matrix with that value along
    %                      the diagonal is returned at each time.
    %                      If R is a vector of same dimension as the 
    %                      observations, a matrix with that vector on the 
    %                      diagonal is returned at each time.
    %                      If R is a matrix, that matrix is returned at
    %                      each time. 
    %                      Otherwise, R is the provided function. 
    %     EX0           -- initial state mean vector (dim_x x 1). 
    %                      Default is zeros(dim_x,1). 
    %     CovX0         -- initial state covariance matrix (dim_x x dim_x).
    %                      Default is diagonal matrix with variance 0.1 
    %                      for bound parameters and variance 0.5 otherwise.
    %     alpha0        -- vector (dim_x x 1) of initial state
    %                      probabilities. Default is
    %                      ones(num_combos,1)/num_combos.
    %
    % 
    % Modified: 
    % Created: 20190601 -- Patrick Stokes
    %
    
properties (SetAccess = private)
        numCombos = [];
        transMatr = [];
        
        offIdxs = [];
        boundIdxs = [];
        
        dimX = [];
        dimY = []
    end
    
    properties (Access = public)
        isDynamic = [];
        peakCombos = [];
        probTrans = [];
        
        f = [];
        Phi = [];
        fIntegr = [];
        fIntegrParams = [];
        Q = [];
        omega = [];
        R = [];
        
        EX0 = [];
        CovX0 = [];
        alpha0
        
    end
    
    methods (Access = public)
        %-----------------------------------
        %       CONSTRUCTOR
        %-----------------------------------
        function obj = StateSpaceMultiPeak(peakModels, isDynamic, peakCombos, probTrans, f, Phi, fIntegr, fIntegrParams, Q, omega, R, EX0, CovX0, alpha0)
            %Constructor for the StateSpaceMultiPeakModel class
            %
            %   obj = StateSpaceMultiPeak(...)
            %
            
            %*************************
            % Handle variable inputs *
            %*************************
            if nargin < 14
                alpha0 = [];
            end
            if nargin < 13
                CovX0 = [];
            end
            if nargin < 12
                EX0 = [];
            end
            if nargin < 11
                R = [];
            end
            if nargin < 10
                omega = [];
            end
            if nargin < 9
                Q = [];
            end
            if nargin < 8
                fIntegrParams = [];
            end
            if nargin < 7
                fIntegr = [];
            end
            if nargin < 6
                Phi = [];
            end
            if nargin < 5
                f = [];
            end
            if nargin < 4
                probTrans = [];
            end
            if nargin < 3
                peakCombos = [];
            end
            if nargin < 2
                isDynamic = [];
            end
            if nargin < 1
                peakModels = [];
            end
            
            %***************************************
            % Construct StateSpaceMultiPeak object *
            %***************************************
            
            % Peaks
            obj@MultiPeakModel(peakModels);
            
            % On/Off-peaks and combinations
            if isempty(isDynamic)
                obj.isDynamic = ones(obj.numPeaks,1);
            else
                obj.isDynamic = isDynamic;
            end
            if isempty(peakCombos)
                obj.peakCombos = makePeakCombos(obj.isDynamic);
            else
                obj.peakCombos = peakCombos;
            end
            obj.numCombos = size(obj.peakCombos,1);
            
            % Pre-compute the indices of peak parameters
            obj.offIdxs = cell(obj.numCombos,1);
            for jj = 1:obj.numCombos
                curr_combo = obj.peakCombos(jj,:);
                obj.offIdxs{jj} = obj.getPeakIdxs(~logical(curr_combo));
            end

            % Combo transition matrix
            if isempty(probTrans)
                obj.probTrans = [0.2 0.2 0.9];
            else
                obj.probTrans = probTrans;
            end
            obj.transMatr = makeComboTransitionMatr(obj.peakCombos,obj.probTrans);
            
            
            obj.dimX = obj.totalPeakParams;
            
            % Drift function
            if ~isempty(f)
                obj.f = f;
            end
            
            % State-transition matrix, integral of diffusion function
            if isempty(Phi)
                obj.Phi = @(t2,t1)(eye(obj.dimX));
            else
                obj.Phi = Phi;
            end
            
            % Integral of drift function and its parameters
            if isempty(fIntegr)
                obj.fIntegr = @(x,t1,t2,params)(zeros(obj.dimX,1));
            else
                obj.fIntegr = fIntegr;
            end
            if isempty(fIntegrParams)
                obj.fIntegrParams = {};
            else
                obj.fIntegrParams = fIntegrParams;
            end
            
            % State noise covariance function
            if isempty(Q)
                obj.Q = @(t)(0.1*eye(obj.dimX));
            elseif isa(Q,'function_handle')
                obj.Q = Q;
            elseif size(Q,1)==1 && size(Q,2)==1 && size(Q,3)==1
                obj.Q = @(t)(Q*eye(obj.dimX));
            elseif ((size(Q,1)==obj.dimX && size(Q,2)==1) || (size(Q,2)==obj.dimX && size(Q,1)==1)) && size(Q,3)==1
                obj.Q = @(t)(diag(Q));
            elseif size(Q,1)==obj.dimX && size(Q,2)==obj.dimX && size(Q,3)==1
                obj.Q = @(t)(Q);
            else
                obj.Q = Q;
            end
            
            % Observation frequencies and dimension
            if isempty(omega)
                obj.omega = linspace(0,100,1000);
            else
                obj.omega = omega;
            end
            obj.dimY = length(obj.omega);
            
            % Observation noise covariance function
            if isempty(R) 
                obj.R = @(t)(eye(obj.dimY));
            elseif isa(R,'function_handle')
                obj.R = R;
            elseif size(R,1)==1 && size(R,2)==1 && size(R,3)==1
                obj.R = @(t)(R*eye(obj.dimY));
            elseif ((size(R,1)==obj.dimY && size(R,2)==1) || (size(R,2)==obj.dimY && size(R,1)==1)) && size(R,3)==1
                obj.R = @(t)(diag(R));
            elseif size(R,1)==obj.dimY && size(R,2)==obj.dimY && size(R,3)==1
                obj.R = @(t)(R);
            else
                obj.R = R;
            end

            % Indexes of bound parameters
            obj.boundIdxs = strcmpi({obj.getParamClasses},'sigmoid');
            
            % Initial state mean
            if isempty(EX0)
                obj.EX0 = zeros(obj.dimX,1);
            else
                obj.EX0 = EX0;
            end
            
            % Initial state covariance
            if isempty(CovX0)
                obj.CovX0 = ones(obj.dimX,1)*.5;
                obj.CovX0(obj.boundIdxs) = .1;
                obj.CovX0 = diag(obj.CovX0);
            else
                obj.CovX0 = CovX0;
            end
            
            % Initial combo probabilities
            if isempty(alpha0)
                obj.alpha0 = ones(obj.numCombos,1)/obj.numCombos;
            else
                obj.alpha0 = alpha0;
            end
            
        end
        
        function [Y,X,alpha,is_on,t,W,V,comps] = simulate(obj,N,f_limit_walk)
            % Function to simulate a StateSpaceMultiPeak object.
            %
            %*** NOTE THIS FUNCTION IS NOT GENERALIZED IT ASSUMES IDENTITY 
            %    TRANSITION, CONSTANT NOISE COVARIANCES, AND TIMES GIVEN 
            %    BY THE INTEGER TIME STEPS.
            %
            % INPUTS:
            %   N            -- number of time steps in simulation 
            %   f_limit_walk -- indicator of whether to hard-bound 
            %                   the states at +/-3. Default 0.
            %
            % OUTPUTS:
            %   Y     -- matrix (dim_y x N) observations
            %   X     -- matrix (dim_x x N+1) of state values 
            %   alpha -- matrix (num_combos x N+1) of combo probabilities
            %   is_on -- indicator matrix (num_peaks x N+1) of which peaks
            %            are on at each time 
            %   t     -- vector of times, 0:N
            %   W     -- matrix (dim_x x N) of state input noises
            %   V     -- matrix (dim_y x N) of observation noises
            %   comps -- array (num_peaks x dim_y x N) of individual peak
            %            components
            %
            
            % Handle variable inputs
            if nargin < 3
                f_limit_walk = [];
            end
            if nargin < 2
                N = [];
            end
            
            % Set defaults
            if isempty(f_limit_walk)
                f_limit_walk = 0;
            end
            if isempty(N)
                N = 100;
            end
            
            % Storage of output variables
            t = 0:N;
            X = zeros(obj.totalPeakParams,N+1);
            alpha = zeros(obj.numCombos,N+1);
            is_on = zeros(obj.numPeaks,N+1);
            Y = zeros(length(obj.omega),N);
            comps = zeros(obj.numPeaks,length(obj.omega),N);
            
            % Initial state
            curr_combo = find(rand <= cumsum(obj.alpha0),1,'first');
            alpha(curr_combo,1) = 1;
            sqrtCovX0 = chol(obj.CovX0);
            x0 = obj.EX0' + randn(1,obj.dimX)*sqrtCovX0;
            X(:,1) = x0';
            
            % Draw state and observation noises
            sqrtQ = chol(obj.Q(1));
            W = randn(N,obj.dimX)*sqrtQ;
            W = W';
            sqrtR = chol(obj.R(1));
            V = randn(N,obj.dimY)*sqrtR;
            V = V';
            
            % Find bounded parameters
            is_exp = strcmpi(obj.getParamClasses,'exponential');
            is_sig = strcmpi(obj.getParamClasses,'sigmoid');

            % Evolve at each time step
            for tt = 2:N+1
                % Update On/Off-combination
                tmp_alpha = obj.transMatr * alpha(:,tt-1);
                curr_combo = find(rand <= cumsum(tmp_alpha),1,'first');
                alpha(curr_combo,tt) = 1;
                % Update state
                X(:,tt) = X(:,tt-1) + W(:,tt-1);
                if f_limit_walk
                    X(is_sig,tt) = max(-3,X(is_sig,tt));
                    X(is_sig,tt) = min(3,X(is_sig,tt));
                    W(is_sig,tt-1) = X(is_sig,tt) - X(is_sig,tt-1);
                    X(is_exp,tt) = min(3,X(is_exp,tt));
                    W(is_exp,tt-1) = X(is_exp,tt) - X(is_exp,tt-1);
                end
                % Form observations
                is_on(:,tt) = obj.peakCombos(curr_combo,:)';
                [~, comps(:,:,tt-1)] = obj.getH(X(:,tt),obj.omega,is_on(:,tt));
                Y(:,tt-1) = comps(:,:,tt-1)'*is_on(:,tt) + V(:,tt-1);
            end
        end
    end
end

