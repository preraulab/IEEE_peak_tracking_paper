classdef PeakModel < handle
    %PeakModel Parametric model of a single spectral peak
    %
    %   obj = PeakModel(modelFunClass, paramLinkFuns, paramLinkFunParams, modelOptions, peakName)
    %
    %   Input:
    %       modelFunClass: string, model function class: 'gaussian', 'box',
    %           'harmonic', 'skewgamma', 'expdecay', 'expdoubledecay',
    %           'smoothdoubledecay', 'linear'
    %       paramLinkFuns: 1 x numPeakParams cell array of strings of link
    %           function classes: 'sigmoid' (default), 'exp', 'identity', 'poly'
    %       paramLinkFunParams: 1 x numPeakParams cell array of
    %           strings of link function parameters
    %       modelOptions: model specific options
    %       peakName: string, descriptive name of peak
    %
    %
    %   PEAK FUNCTIONS:
    %        Gaussian ('gaussian')
    %             Function:
    %                 y(x) = A * exp(-(w-F)^2 / (2*B))
    %             modelParams:
    %                 F - central frequency
    %                 A - maximum amplitude
    %                 B - bandwidth (variance)
    %
    %         Box Function ('box');
    %             Function:
    %                 y(x) = A * exp(-(w-F)^p / (2*B))
    %             modelParams:
    %                 F - central frequency
    %                 A - maximum amplitude
    %                 B - bandwidth (variance)
    %             modelOptions:
    %                 p - polynomial order
    %
    %         Gaussian with Harmonics ('harmonic');
    %             Function:
    %                 y(x) = SUM h = 0 to H_peak: (Beta^h) * A * exp(-(w - (h+1)*F)^2 / (2*B));
    %             modelParams:
    %                 F - central frequency
    %                 A - maximum amplitude
    %                 B - bandwidth (variance)
    %                 Beta - harmonic ampiltude multiplier
    %             modelOptions:
    %                 H_peak - number of harmonics
    %
    %         Gamma with Skewness ('skewgamma');
    %             Function:
    %                 A gamma function with fixed mode, variance, and skew
    %             modelParams:
    %                 F - central frequency (mode)
    %                 A - maximum amplitude
    %                 B - bandwidth (variance)
    %                 S - skewness
    %
    %         Exponential Decay ('expdecay');
    %             Function:
    %                  y(x) = BA.*(1-BR).^x +BO;
    %             modelParams:
    %                 BA - maximum amplitude
    %                 BR - decay rate
    %                 BO - offset/baseline
    %     
    %     *****************************************
    %     SECOND DERIVATIVES NOT IMPLEMENTED YET
    %     *****************************************
    %         Double Exponential Decay ('expdecay');
    %             Function:
    %                  y(x) = BA1.*(1-BR1).^x +  BA2.*(1-BR2).^x +BO;
    %             modelParams:
    %                 BA1 - maximum amplitude of first exponential
    %                 BR1 - decay rate of first exponential
    %                 BA2 - maximum amplitude of second exponential
    %                 BR2 - decay rate of second exponential
    %                 BO - offset/baseline
    %
    %         Piecewise Exponential with Smooth Transition ('smoothdoubledecay');
    %             Function:
    %                  see: Zelterman et. al, "Piecewise exponential survival curves with
    %                  smooth transitions", Mathematical Biosciences, 1994
    %             modelParams:
    %                 L1mL2 - difference between the rates of the two
    %                 exponentials
    %                 L2 - The rate of the second exponential
    %                 K - maximum amplitude
    %                 O - offset/baseline
    %
    %   Copyright 2019 Prerau Laboratory
    %
    %   Last modified 5/24/2019 by Michael Prerau
    
    properties (Access = public)
        modelFunClass; %Function class for peak model
        
        modelParams = PeakObjParam; %List of state parameters
        
        paramNames; %State variable name
        paramDescriptions; %State variable description
        
        modelOptions; %Additional variables per model
        
        modelFun; %Handle to model function
        modelDerivFuns; %Cell array of partial derivative function handles
        model2ndDerivFuns; %Cell array of partial derivative function handles
        
        peakName; %String name for peak
    end
    
    properties (SetAccess = private)
        numPeakParams; %Number of peak parameters (cannot set)
    end
    
    methods (Access = public)
        %-----------------------------------
        %       CONSTRUCTOR
        %-----------------------------------
        function obj = PeakModel(modelFunClass, paramLinkFuns, paramLinkFunParams, modelOptions, peakName)
            %Constructor for the PeakModel class
            %
            %   obj = PeakModel(modelFunClass, paramLinkFuns, paramLinkFunParams, modelOptions, peakName)
            %
            %   Input:
            %       modelFunClass: string, model function class: 'gaussian', 'box',
            %           'harmonic', 'skewgamma', 'expdecay', 'expdoubledecay',
            %           'smoothdoubledecay', 'linear'
            %       paramLinkFuns: 1 x numPeakParams cell array of strings of link
            %           function classes: 'sigmoid' (default), 'exp', 'identity', 'poly'
            %       paramLinkFunParams: 1 x numPeakParams cell array of
            %           strings of link function parameters
            %       modelOptions: model specific options
            %       peakName: string, descriptive name of peak
            
            if nargin == 0
                modelFunClass = 'Gaussian';
            end
            
            if nargin < 2
                paramLinkFuns = [];
            end
            
            if nargin < 3
                paramLinkFunParams = [];
            end
            
            if nargin < 4
                modelOptions = [];
            end
            
            if nargin < 5
                peakName = [];
            end
            
            obj.modelOptions = modelOptions;
            
            switch lower(modelFunClass)
                %Sigmoid link function
                case {'gaussian','normal'}
                    obj.modelFunClass = 'Gaussian';
                    
                    %Define model function
                    obj.modelFun = @obj.gaussian;
                    obj.paramDescriptions = {'Peak Frequency', 'Amplitude', 'Bandwidth'};
                    
                    %Define model partial derivative functions
                    obj.modelDerivFuns = {@obj.gaussian_dH_dFbar, @obj.gaussian_dH_dAbar, @obj.gaussian_dH_dBbar};
                    
                    %Define model second partial derivative functions
                    obj.model2ndDerivFuns = ...
                        {@obj.gaussian_d2H_dFbar2,      @obj.gaussian_d2H_dFbar_dAbar,  @obj.gaussian_d2H_dFbar_dBbar; ...
                        @obj.gaussian_d2H_dFbar_dAbar,  @obj.gaussian_d2H_dAbar2,       @obj.gaussian_d2H_dAbar_dBbar; ...
                        @obj.gaussian_d2H_dFbar_dBbar,  @obj.gaussian_d2H_dAbar_dBbar,  @obj.gaussian_d2H_dBbar2};
                    
                    %Set default model link functions
                    if isempty(paramLinkFuns)
                        paramLinkFuns = {'Sigmoid', 'Sigmoid', 'Sigmoid'};
                    end
                    
                    %Set default model link function parameters
                    if isempty(paramLinkFunParams)
                        paramLinkFunParams = {[0 1], [0 1], [0 1]};
                    end
                    
                    %Create the parameters
                    obj.modelParams(1) =  PeakObjParam('F', paramLinkFuns{1}, paramLinkFunParams{1});
                    obj.modelParams(2) =  PeakObjParam('A', paramLinkFuns{2}, paramLinkFunParams{2});
                    obj.modelParams(3) =  PeakObjParam('B', paramLinkFuns{3}, paramLinkFunParams{3});
                    
                    
                case {'gaussian harmonic','harmonic', 'gaussian_harmonic', 'normal harmonic','normal_harmonic'}
                    obj.modelFunClass = 'Gaussian Harmonic';
                    
                    %Define model function
                    obj.modelFun = @obj.gaussian_harmonic;
                    obj.paramDescriptions = {'Peak Frequency', 'Amplitude', 'Bandwidth', 'Harmonic Ratio'};
                    
                    %Set two harmonics by default
                    if isempty(obj.modelOptions)
                        obj.modelOptions = 2;
                    end
                    
                    %Define model partial derivative functions
                    obj.modelDerivFuns = {@obj.gaussian_harmonic_dH_dFbar, @obj.gaussian_harmonic_dH_dAbar, @obj.gaussian_harmonic_dH_dBbar , @obj.gaussian_harmonic_dH_dBetabar};
                    
                    %Define model second partial derivative functions
                    obj.model2ndDerivFuns = ...
                        {@obj.gaussian_harmonic_d2H_dFbar2,         @obj.gaussian_harmonic_d2H_dFbar_dAbar,     @obj.gaussian_harmonic_d2H_dFbar_dBbar,     @obj.gaussian_harmonic_d2H_dFbar_dBetabar; ...
                        @obj.gaussian_harmonic_d2H_dFbar_dAbar,     @obj.gaussian_harmonic_d2H_dAbar2,          @obj.gaussian_harmonic_d2H_dAbar_dBbar,     @obj.gaussian_harmonic_d2H_dAbar_dBetabar; ...
                        @obj.gaussian_harmonic_d2H_dFbar_dBbar,     @obj.gaussian_harmonic_d2H_dAbar_dBbar,     @obj.gaussian_harmonic_d2H_dBbar2,          @obj.gaussian_harmonic_d2H_dBbar_dBetabar; ...
                        @obj.gaussian_harmonic_d2H_dFbar_dBetabar,  @obj.gaussian_harmonic_d2H_dAbar_dBetabar,  @obj.gaussian_harmonic_d2H_dBbar_dBetabar,  @obj.gaussian_harmonic_d2H_dBetabar2};
                    
                    %Set default model link functions
                    if isempty(paramLinkFuns)
                        paramLinkFuns = {'Sigmoid', 'Sigmoid', 'Sigmoid', 'Sigmoid'};
                    end
                    
                    %Set default model link function parameters
                    if isempty(paramLinkFunParams)
                        paramLinkFunParams = {[0 1], [0 1], [0 1], [0 1]};
                    end
                    
                    %Create the parameters
                    obj.modelParams(1) =  PeakObjParam('F',paramLinkFuns{1},paramLinkFunParams{1});
                    obj.modelParams(2) =  PeakObjParam('A',paramLinkFuns{2},paramLinkFunParams{2});
                    obj.modelParams(3) =  PeakObjParam('B',paramLinkFuns{3},paramLinkFunParams{3});
                    obj.modelParams(4) =  PeakObjParam('Beta',paramLinkFuns{4},paramLinkFunParams{4});
                    
                case {'box', 'smoothed box'}
                    obj.modelFunClass = 'Box';
                    
                    %Define model function
                    obj.modelFun = @obj.boxfun;
                    obj.paramDescriptions = {'Peak Frequency', 'Amplitude', 'Bandwidth'};
                    
                    %Set default for p
                    if isempty(obj.modelOptions)
                        obj.modelOptions = 16;
                    end
                    
                    %Define model partial derivative functions
                    obj.modelDerivFuns = {@obj.boxfun_dH_dFbar, @obj.boxfun_dH_dAbar, @obj.boxfun_dH_dBbar};
                    
                    %Define model second partial derivative functions
                    obj.model2ndDerivFuns = ...
                        {@obj.box_d2H_dFbar2,      @obj.box_d2H_dFbar_dAbar,  @obj.box_d2H_dFbar_dBbar; ...
                        @obj.box_d2H_dFbar_dAbar,  @obj.box_d2H_dAbar2,       @obj.box_d2H_dAbar_dBbar; ...
                        @obj.box_d2H_dFbar_dBbar,  @obj.box_d2H_dAbar_dBbar,  @obj.box_d2H_dBbar2};
                    
                    %Set default model link functions
                    if isempty(paramLinkFuns)
                        paramLinkFuns = {'Sigmoid', 'Sigmoid', 'Sigmoid'};
                    end
                    
                    %Set default model link function parameters
                    if isempty(paramLinkFunParams)
                        paramLinkFunParams = {[0 1], [0 1], [0 1]};
                    end
                    
                    %Create the parameters
                    obj.modelParams(1) =  PeakObjParam('F',paramLinkFuns{1},paramLinkFunParams{1});
                    obj.modelParams(2) =  PeakObjParam('A',paramLinkFuns{2},paramLinkFunParams{2});
                    obj.modelParams(3) =  PeakObjParam('B',paramLinkFuns{3},paramLinkFunParams{3});
                    
                case {'gamma', 'skewgamma','skew_gamma','skew gamma', 'skewed gamma', 'skewed_gamma'}
                    obj.modelFunClass = 'Skewed Gamma';
                    
                    %Define model function
                    obj.modelFun = @obj.skewgamma;
                    obj.paramDescriptions = {'Peak Frequency', 'Amplitude', 'Bandwidth', 'Skewness'};
                    
                    %Define model partial derivative functions
                    obj.modelDerivFuns = {@obj.skewgamma_dH_dFbar, @obj.skewgamma_dH_dAbar, @obj.skewgamma_dH_dBbar , @obj.skewgamma_dH_dSbar};
                    
                    %Define model second partial derivative functions
                    obj.model2ndDerivFuns = ...
                        {@obj.skewgamma_d2H_dFbar2,         @obj.skewgamma_d2H_dFbar_dAbar,     @obj.skewgamma_d2H_dFbar_dBbar,     @obj.skewgamma_d2H_dFbar_dSbar; ...
                        @obj.skewgamma_d2H_dFbar_dAbar,     @obj.skewgamma_d2H_dAbar2,          @obj.skewgamma_d2H_dAbar_dBbar,     @obj.skewgamma_d2H_dAbar_dSbar; ...
                        @obj.skewgamma_d2H_dFbar_dBbar,     @obj.skewgamma_d2H_dAbar_dBbar,     @obj.skewgamma_d2H_dBbar2,          @obj.skewgamma_d2H_dBbar_dSbar; ...
                        @obj.skewgamma_d2H_dFbar_dSbar,     @obj.skewgamma_d2H_dAbar_dSbar,     @obj.skewgamma_d2H_dBbar_dSbar,     @obj.skewgamma_d2H_dSbar2};
                    
                    %Set default model link functions
                    if isempty(paramLinkFuns)
                        paramLinkFuns = {'Sigmoid', 'Sigmoid', 'Sigmoid', 'Sigmoid'};
                    end
                    
                    %Set default model link function parameters
                    if isempty(paramLinkFunParams)
                        paramLinkFunParams = {[0 1], [0 1], [0 1], [.5 1]};
                    end
                    
                    %Create the parameters
                    obj.modelParams(1) =  PeakObjParam('F',paramLinkFuns{1},paramLinkFunParams{1});
                    obj.modelParams(2) =  PeakObjParam('A',paramLinkFuns{2},paramLinkFunParams{2});
                    obj.modelParams(3) =  PeakObjParam('B',paramLinkFuns{3},paramLinkFunParams{3});
                    obj.modelParams(4) =  PeakObjParam('S',paramLinkFuns{4},paramLinkFunParams{4});
                    
                case {'expdecay', 'decay','exponential decay', 'exponential_decay'}
                    obj.modelFunClass = 'Exponential Decay';
                    
                    %Define model function
                    obj.modelFun = @obj.expdecay;
                    obj.paramDescriptions = {'Y-Intercept', 'Rate', 'Offset'};
                    
                    %Define model partial derivative functions
                    obj.modelDerivFuns = {@obj.expdecay_dHdBA_bar, @obj.expdecay_dHdBR_bar, @obj.expdecay_dHdBO_bar};
                    
                                        %Define model second partial derivative functions
                    obj.model2ndDerivFuns = ...
                        {@obj.expdecay_d2H_dBAbar2,      @obj.expdecay_d2H_dBAbar_dBRbar,  @obj.expdecay_d2H_dBAbar_dBObar; ...
                        @obj.expdecay_d2H_dBAbar_dBRbar,  @obj.expdecay_d2H_dBRbar2,       @obj.expdecay_d2H_dBRbar_dBObar; ...
                        @obj.expdecay_d2H_dBAbar_dBObar,  @obj.expdecay_d2H_dBRbar_dBObar,  @obj.expdecay_d2H_dBObar2};
                    
                    %Set default model link functions
                    if isempty(paramLinkFuns)
                        paramLinkFuns = {'Sigmoid', 'Sigmoid', 'Identity'};
                    end
                    
                    %Set default model link function parameters
                    if isempty(paramLinkFunParams)
                        paramLinkFunParams = {[0 10], [0 1], [ ]};
                    end
                    
                    %Create the parameters
                    obj.modelParams(1) =  PeakObjParam('BA',paramLinkFuns{1},paramLinkFunParams{1});
                    obj.modelParams(2) =  PeakObjParam('BR',paramLinkFuns{2},paramLinkFunParams{2});
                    obj.modelParams(3) =  PeakObjParam('BO',paramLinkFuns{3},paramLinkFunParams{3});
                    
                case {'doubleexpdecay', 'doubledecay', 'double decay', 'double_decay','double exponential decay', 'double_exponential_decay'}
                    obj.modelFunClass = 'Exponential Decay';
                    
                    %Define model function
                    obj.modelFun = @obj.doubledecay;
                    obj.paramDescriptions = {'Y-Intercept 1', 'Rate 1', 'Y-Intercept 2', 'Rate 2', 'Offset'};
                    
                    %Define model partial derivative functions
                    obj.modelDerivFuns = {@obj.doubledecay_dHdBA1_bar, @obj.doubledecay_dHdBR1_bar,@obj.doubledecay_dHdBA2_bar, @obj.doubledecay_dHdBR2_bar,@obj.doubledecay_dHdBO_bar};
                    
                    %Set default model link functions
                    if isempty(paramLinkFuns)
                        paramLinkFuns = {'Sigmoid', 'Sigmoid', 'Sigmoid', 'Sigmoid', 'Identity'};
                    end
                    
                    %Set default model link function parameters
                    if isempty(paramLinkFunParams)
                        paramLinkFunParams = {[0 10], [0 1], [0 10], [0 1], [ ]};
                    end
                    
                    %Create the parameters
                    obj.modelParams(1) =  PeakObjParam('BA1',paramLinkFuns{1},paramLinkFunParams{1});
                    obj.modelParams(2) =  PeakObjParam('BR1',paramLinkFuns{2},paramLinkFunParams{2});
                    obj.modelParams(3) =  PeakObjParam('BA2',paramLinkFuns{3},paramLinkFunParams{3});
                    obj.modelParams(4) =  PeakObjParam('BR2',paramLinkFuns{4},paramLinkFunParams{4});
                    obj.modelParams(5) =  PeakObjParam('BO',paramLinkFuns{5},paramLinkFunParams{5});
                case {'gamma double decay', 'smooth double decay', 'gamma_double_decay', 'gamma_double', 'smoothdoubledecay','gammadoubledecay'}
                    obj.modelFunClass = 'Double Exponential Decay with Gamma Transition';
                    
                    %Define model function
                    obj.modelFun = @obj.smoothdoubledecay;
                    obj.paramDescriptions = {'L1 minus L2', 'L2', 'K', 'Offset'};
                    
                    %Define model partial derivative functions
                    obj.modelDerivFuns = {@obj.smoothdoubledecay_dHdL1mL2t_bar, @obj.smoothdoubledecay_dHdL2_bar, ...
                        @obj.smoothdoubledecay_dHdK_bar, @obj.smoothdoubledecay_dHdO_bar};
                    
                    %Set default model link functions
                    if isempty(paramLinkFuns)
                        paramLinkFuns = {'Sigmoid', 'Sigmoid', 'Sigmoid', 'Identity'};
                    end
                    
                    %Set default model link function parameters
                    if isempty(paramLinkFunParams)
                        paramLinkFunParams = {[0 10], [0 1], [0 10], [0 1], [ ]};
                    end
                    
                    %Create the parameters
                    obj.modelParams(1) =  PeakObjParam('L1mL2',paramLinkFuns{1},paramLinkFunParams{1});
                    obj.modelParams(2) =  PeakObjParam('L2',paramLinkFuns{2},paramLinkFunParams{2});
                    obj.modelParams(3) =  PeakObjParam('K',paramLinkFuns{3},paramLinkFunParams{3});
                    obj.modelParams(4) =  PeakObjParam('O',paramLinkFuns{4},paramLinkFunParams{4});
                    
                    %Set default for p
                    if isempty(obj.modelOptions)
                        obj.modelOptions = [1 2];
                    end
                    
                case {'oof'}
                    obj.modelFunClass = 'One-Over-F';
                    
                    %Define model function
                    obj.modelFun = @obj.oof;
                    obj.paramDescriptions = {'Knee', 'Rate', 'Baseline Offset'};
                    
                    %Define model partial derivative functions
                    obj.modelDerivFuns = {@obj.oof_dHdK_bar, @obj.oof_dHdR_bar, @obj.oof_dHdB_bar};
                    
                                        %Define model second partial derivative functions
                    obj.model2ndDerivFuns = ...
                        {@obj.oof_d2H_dKbar2,      @obj.oof_d2H_dKbar_dRbar,  @obj.oof_d2H_dKbar_dBbar; ...
                        @obj.oof_d2H_dKbar_dRbar,  @obj.oof_d2H_dRbar2,       @obj.oof_d2H_dRbar_dBbar; ...
                        @obj.oof_d2H_dKbar_dBbar,  @obj.oof_d2H_dRbar_dBbar,  @obj.oof_d2H_dBbar2};
                    
                    %Set default model link functions
                    if isempty(paramLinkFuns)
                        paramLinkFuns = {'Sigmoid', 'Sigmoid', 'Identity'};
                    end
                    
                    %Set default model link function parameters
                    if isempty(paramLinkFunParams)
                        paramLinkFunParams = {[0 10], [0 1], [ ]};
                    end
                    
                    %Create the parameters
                    obj.modelParams(1) =  PeakObjParam('K',paramLinkFuns{1},paramLinkFunParams{1});
                    obj.modelParams(2) =  PeakObjParam('R',paramLinkFuns{2},paramLinkFunParams{2});
                    obj.modelParams(3) =  PeakObjParam('B',paramLinkFuns{3},paramLinkFunParams{3});
                    
                case {'linear'}
                    obj.modelFunClass = 'Linear';
                    
                    %Define model function
                    obj.modelFun = @obj.linear;
                    obj.paramDescriptions = {'Slope', 'Y-Intercept'};
                    
                    %Define model partial derivative functions
                    obj.modelDerivFuns = {@obj.linear_dH_dMbar, @obj.linear_dH_dBbar};
                    
                    %Define model second partial derivative functions
                    obj.model2ndDerivFuns = ...
                        {@obj.linear_dH_dMbar2,       @obj.linear_dH_dMbar_dBbar; ...
                        @obj.linear_dH_dMbar_dBbar,   @obj.linear_dH_dBbar2};
                    
                    %Set default model link functions
                    if isempty(paramLinkFuns)
                        paramLinkFuns = {'Sigmoid', 'Identity'};
                    end
                    
                    %Set default model link function parameters
                    if isempty(paramLinkFunParams)
                        paramLinkFunParams = {[0 10], [ ]};
                    end
                    
                    %Create the parameters
                    obj.modelParams(1) =  PeakObjParam('M',paramLinkFuns{1},paramLinkFunParams{1});
                    obj.modelParams(2) =  PeakObjParam('B',paramLinkFuns{2},paramLinkFunParams{2});
                otherwise
                    error('Invalid model class');
            end
            
            %Number of parameters
            obj.numPeakParams = length(obj.modelParams);
            obj.paramNames = {obj.modelParams.paramName};
            obj.peakName = peakName;
        end
        
        %-----------------------------------
        % PLOT PEAK FUNCTION
        %-----------------------------------
        function plot(obj, X, w)
            %Plot function for PeakModel class
            %
            %  plot(obj, X, w, <plot function arguments>)
            %
            %   Input:
            %       X: 1 x totalPeakParams double, state vector (default: zeros)
            %       w: 1 x W double, frequency vector
            
            if nargin < 3 || isempty(X)
                X=zeros(1, obj.numPeakParams);
            end
            
            if nargin < 2 || isempty(w)
                %Figure out the range of frequencies the covers the peak
                freq_idx = find(strcmpi({obj.modelParams.paramName},'F'));
                bw_idx = find(strcmpi({obj.modelParams.paramName},'B'));
                
                %Set the range with a default of 0-100Hz for non-peak functions
                if ~isempty(freq_idx) && ~isempty(bw_idx)
                    F = obj.modelParams(freq_idx).linkFun(X(freq_idx));
                    B = obj.modelParams(bw_idx).linkFun(X(freq_idx));
                    
                    if strcmpi(obj.modelFunClass, 'gaussian harmonic')
                        w=linspace(F - 4*B, F*obj.modelOptions + 4*B, 10000);
                    else
                        w=linspace(F - 4*B, F + 4*B, 10000);
                    end
                else
                    w=linspace(0, 100, 10000);
                end
            end
            
            H_peak = obj.modelFun(w,X);
            plot(w, H_peak);
        end
        
        %-----------------------------------
        % NUMERICAL DERIVATIVES
        %-----------------------------------
        
        function testDerivs(obj, X, dX)
            %Test to see if the derivatives are correct against numerical approximation
            %
            %  testDerivs(dX)
            %
            %   Input:
            %       dX: double, size of state difference (default: 1e-3)
            
            %Set default for dX
            if nargin<3
                dX =[];
            end
            if nargin<2
                X = [];
            end
            
            if isempty(dX)
                dX=1e-3;
            end
            if isempty(X)
                X=zeros(obj.numPeakParams, 1);
            end
            
            %Figure out the range of frequencies the covers the peak
            freq_idx = find(strcmpi({obj.modelParams.paramName},'F'));
            bw_idx = find(strcmpi({obj.modelParams.paramName},'B'));
            
            %Set the range with a default of 0-100Hz for non-peak functions
            if ~isempty(freq_idx) && ~isempty(bw_idx)
                F = obj.modelParams(freq_idx).linkFun(0);
                B = obj.modelParams(bw_idx).linkFun(0);
                
                if strcmpi(obj.modelFunClass, 'gaussian harmonic')
                    w=linspace(F - 4*B, F*obj.modelOptions + 4*B, 10000);
                else
                    w=linspace(F - 4*B, F + 4*B, 10000);
                end
            else
                w=linspace(0, 100, 10000);
            end
            
            %Loop over each parameter an plot comparison with numerical difference
            for ii = 1:obj.numPeakParams
                %Create X and adjusted X
                Xdp = X;
                Xdp(ii) = Xdp(ii)+dX;
                
                %Compute both H_peak functions and take difference
                H_peak = obj.modelFun(w,X);
                H_peak_2 = obj.modelFun(w,Xdp);
                
                emp_deriv = (H_peak_2-H_peak)/dX;
                
                analytic_deriv = obj.modelDerivFuns{ii}(w,X, H_peak);
                
                %Plot results
                figure
                hold all
                plot(w, H_peak, 'linewidth',2)
                plot(w, analytic_deriv, 'linewidth',2);
                plot(w, emp_deriv, '--', 'linewidth',2);
                
                legend('Function', 'Analytic', 'Numerical');
                
                title(obj.modelParams(ii).paramName);
            end
        end
        
        function test2ndDerivs(obj,X,dX)
            %Test to see if the derivatives are correct against numerical approximation
            %
            %  test2ndDerivs(dX)
            %
            %   Input:
            %       dX: double, size of state difference (default: 1e-3)
            
            %Set default for dX
            if nargin<3
                dX = [];
            end
            if nargin<2
                X = [];
            end
            
            if isempty(dX)
                dX=1e-5;
            end
            
            if isempty(X)
                X=zeros(obj.numPeakParams, 1);
            end
            
            %Figure out the range of frequencies the covers the peak
            freq_idx = find(strcmpi({obj.modelParams.paramName},'F'));
            bw_idx = find(strcmpi({obj.modelParams.paramName},'B'));
            
            %Set the range with a default of 0-100Hz for non-peak functions
            if ~isempty(freq_idx) && ~isempty(bw_idx)
                F = obj.modelParams(freq_idx).linkFun(0);
                B = obj.modelParams(bw_idx).linkFun(0);
                
                if strcmpi(obj.modelFunClass, 'gaussian harmonic')
                    w=linspace(F - 4*B, F*obj.modelOptions + 4*B, 10000);
                else
                    w=linspace(F - 4*B, F + 4*B, 10000);
                end
            else
                w=linspace(0, 100, 10000);
            end
            
            
            H_peak = obj.modelFun(w,X);
            
            %Loop over each parameter an plot comparison with numerical difference
            for zz1 = 1:obj.numPeakParams
                for zz2 = zz1:obj.numPeakParams
                    
                    analytic_2nd_deriv = obj.model2ndDerivFuns{zz1,zz2}(w,X, H_peak);
                    
                    dx1 = zeros(size(X));
                    dx1(zz1) = dX;
                    dx2 = zeros(size(X));
                    dx2(zz2) = dX;
                    
                    if zz1~=zz2
                        hx1px2p = obj.modelFun(w,X+dx1+dx2);
                        hx1px2m = obj.modelFun(w,X+dx1-dx2);
                        hx1mx2p = obj.modelFun(w,X-dx1+dx2);
                        hx1mx2m = obj.modelFun(w,X-dx1-dx2);
                        emp_2nd_deriv = (hx1px2p - hx1px2m - hx1mx2p + hx1mx2m)/(4*dX*dX);
                    else
                        hx1px2 = obj.modelFun(w,X+dx1);
                        hx1mx2 = obj.modelFun(w,X-dx1);
                        emp_2nd_deriv = (hx1px2 - 2*H_peak + hx1mx2)/(dX*dX);
                    end
                    
                    
                    %Plot results
                    figure
                    hold all
                    plot(w, H_peak, 'linewidth',2)
                    plot(w, analytic_2nd_deriv, 'linewidth',2);
                    plot(w, emp_2nd_deriv, '--', 'linewidth',2);
                    
                    legend('Function', 'Analytic', 'Numerical');
                    
                    title([obj.modelFunClass ': d/d' obj.modelParams(zz1).paramName ' d/d' obj.modelParams(zz2).paramName]);
                    
                end
            end
        end
    end
    
    methods (Access = private)
        %% GAUSSIAN MODEL
        
        %-----------------------------------
        % GAUSSIAN PEAK MODEL
        %-----------------------------------
        function H_peak = gaussian(obj, w, X)
            if size(X,1) ~= 3
                error('Gaussian model requires 3 parameters: F, A, and B');
            end
            
            %Get param values through link functions
            F=obj.modelParams(1).linkFun(X(1,:));
            A=obj.modelParams(2).linkFun(X(2,:));
            B=obj.modelParams(3).linkFun(X(3,:));
            
            H_peak=A.*exp(-(w-F).^2./(2*B));
        end
        
        %-----------------------------------
        % GAUSSIAN PEAK MODEL DERIVAIVES
        %-----------------------------------
        
        function out = gaussian_dH_dFbar(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            
            out = H_peak.*(w-F)/(B).*dF_dFbar;
        end
        
        function out = gaussian_dH_dAbar(obj, ~, X, H_peak)
            %Get param values through link functions
            A = obj.modelParams(2).linkFun(X(2,:));
            
            %Get derivatives of link functions
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            
            out = H_peak./A.*dA_dAbar;
        end
        
        function out = gaussian_dH_dBbar(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            out = H_peak.*(w-F).^2./(2*B.^2).*dB_dBbar;
        end
        
        %-----------------------------------
        % GAUSSIAN PEAK MODEL 2nd DERIVAIVES
        %-----------------------------------
        function out = gaussian_d2H_dFbar2(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            
            Z = (w-F)./B;
            dZ_dF = -1/B;
            out = H_peak.*( (Z.^2 + dZ_dF).*(dF_dFbar.^2) + Z.*d2F_dFbar2);
        end
        
        function out = gaussian_d2H_dFbar_dAbar(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            
            out = H_peak.*((w-F)./B).*(1/A).*dF_dFbar.*dA_dAbar;
        end
        
        function out = gaussian_d2H_dFbar_dBbar(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            
            out = H_peak.*(((w-F).^3)./(2*B.^3)-(w-F)./(B.^2)).*dF_dFbar.*dB_dBbar;
        end
        
        function out = gaussian_d2H_dAbar2(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            
            out = H_peak.*(1/A).*d2A_dAbar2;
        end
        
        function out = gaussian_d2H_dAbar_dBbar(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            
            out = H_peak.*( ((w-F).^2)./(2*B.^2) .* (1/A) ) .* dA_dAbar .* dB_dBbar;
        end
        
        function out = gaussian_d2H_dBbar2(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            
            Z = ((w-F).^2)./(2*B.^2);
            dZ_dB = -((w-F).^2)./(B.^3);
            
            out = H_peak.*((Z.^2 + dZ_dB).*(dB_dBbar.^2) + Z.*d2B_dBbar2);
        end
        
        %% GAUSSIAN HARMONIC MODEL
        
        %-----------------------------------
        % GAUSSIAN HARMONIC PEAK MODEL
        %-----------------------------------
        function H_peak = gaussian_harmonic(obj, w, X)
            if size(X,1) ~= 4
                error('Gaussian harmonic model requires 4 parameters: F, A, B and Beta');
            end
            
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            Beta = obj.modelParams(4).linkFun(X(4,:));
            
            %Get number of harmonics
            N = obj.modelOptions;
            
            %No weighting on the fundamental
            H_peak = zeros(size(w));
            
            for harmonic_num=0:N
                H_peak = H_peak + (Beta.^harmonic_num).*A.*exp(-(w-(harmonic_num+1)*F).^2./(2*B));
            end
        end
        
        %---------------------------------------------
        % GAUSSIAN HARMONIC PEAK MODEL DERIVAIVES
        %---------------------------------------------
        function out=gaussian_harmonic_dH_dFbar(obj, w, X, ~)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            Beta = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            
            %Get number of harmonics
            N = obj.modelOptions;
            
            out= zeros(size(w));
            
            for harmonic_num = 0:N
                Hh = A.* Beta.^harmonic_num .* exp(-(w-(harmonic_num+1).*F).^2./(2*B));
                out = out+Hh.*((harmonic_num+1).*(w-F.*(harmonic_num+1)))./B;
            end
            
            out=out.*dF_dFbar;
        end
        
        function out=gaussian_harmonic_dH_dAbar(obj, ~, X, H_peak)
            %Get param values through link functions
            A = obj.modelParams(2).linkFun(X(2,:));
            
            %Get derivatives of link functions
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            
            out=H_peak./A.*dA_dAbar;
        end
        
        function out=gaussian_harmonic_dH_dBbar(obj, w, X, ~)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            Beta = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get number of harmonics
            N = obj.modelOptions;
            
            %No weighting on the fundamental
            H0=A.*exp(-(w-F).^2./(2*B));
            
            out= H0.*((w-F).^2)/(2*B.^2);
            
            for harmonic_num=1:N
                Hh=A.* Beta.^harmonic_num .* exp(-(w-(harmonic_num+1).*F).^2./(2*B));
                out=out+Hh.*(w-F.*(harmonic_num+1)).^2./(2*B.^2);
            end
            
            out=out.*dB_dBbar;
            
        end
        
        function out=gaussian_harmonic_dH_dBetabar(obj, w, X, ~)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            Beta = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dBetAdBetA_bar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get number of harmonics
            N = obj.modelOptions;
            
            out=0;
            
            for harmonic_num=1:N
                Hh=A.* Beta.^harmonic_num .* exp(-(w-(harmonic_num+1).*F).^2./(2*B));
                out=out+Hh.*(harmonic_num./Beta);
            end
            
            out=out.*dBetAdBetA_bar;
        end
        
        %---------------------------------------------
        % GAUSSIAN HARMONIC PEAK MODEL 2nd DERIVAIVES
        %---------------------------------------------
        function out=gaussian_harmonic_d2H_dFbar2(obj, w, X, ~)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            Beta = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dBeta_dBetabar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2Beta_dBetabar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Get number of harmonics
            N = obj.modelOptions;
            
            H0 = A.*exp(-(w-F).^2./(2*B));
            out = H0.*(((w-F)./B).^2-1/B);
            first_deriv = H0.*((w-F)./B);
            for nn = 1:N
                Hn = A.* Beta.^nn .* exp(-(w-(nn+1).*F).^2./(2*B));
                out = out + Hn.*(((nn+1).*(w-(nn+1).*F)./B).^2-((nn+1)^2)/B);
                first_deriv = first_deriv + Hn.*((nn+1).*(w-(nn+1).*F)./B);
            end
            
            out = out.*(dF_dFbar.^2) + first_deriv.*d2F_dFbar2;
        end
        function out=gaussian_harmonic_d2H_dFbar_dAbar(obj, w, X, ~)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            Beta = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dBeta_dBetabar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2Beta_dBetabar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Get number of harmonics
            N = obj.modelOptions;
            
            H0 = A.*exp(-(w-F).^2./(2*B));
            out = H0.*((w-F)./B).*(1/A);
            for nn = 1:N
                Hn = A.* Beta.^nn .* exp(-(w-(nn+1).*F).^2./(2*B));
                out = out + Hn.*((nn+1).*(w-(nn+1).*F)./B).*(1/A);
            end
            
            out = out.*dF_dFbar.*dA_dAbar;
        end
        
        function out=gaussian_harmonic_d2H_dFbar_dBbar(obj, w, X, ~)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            Beta = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dBeta_dBetabar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2Beta_dBetabar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Get number of harmonics
            N = obj.modelOptions;
            
            H0 = A.*exp(-(w-F).^2./(2*B));
            out = H0.*(((w-F).^3)./(2*B.^3)-(w-F)./(B.^2));
            for nn = 1:N
                Hn = A.* Beta.^nn .* exp(-(w-(nn+1).*F).^2./(2*B));
                out = out + Hn.*(((nn+1)*(w-(nn+1)*F).^3)./(2*B.^3)-((nn+1)*(w-(nn+1)*F))./(B.^2));
            end
            
            out = out.*dF_dFbar.*dB_dBbar;
        end
        
        function out=gaussian_harmonic_d2H_dFbar_dBetabar(obj, w, X, ~)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            Beta = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dBeta_dBetabar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2Beta_dBetabar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Get number of harmonics
            N = obj.modelOptions;
            
            H0 = A.*exp(-(w-F).^2./(2*B));
            out = H0.*( (w-F)./B .* (0/Beta) ); %**Note this is zero.
            for nn = 1:N
                Hn = A.* Beta.^nn .* exp(-(w-(nn+1).*F).^2./(2*B));
                out = out + Hn.*( (nn+1)*(w-(nn+1)*F)./B .* (nn/Beta) );
            end
            
            out = out.*dF_dFbar.*dBeta_dBetabar;
        end
        
        function out=gaussian_harmonic_d2H_dAbar2(obj, w, X, ~)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            Beta = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dBeta_dBetabar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2Beta_dBetabar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Get number of harmonics
            N = obj.modelOptions;
            
            H0 = A.*exp(-(w-F).^2./(2*B));
            out = H0;
            for nn = 1:N
                Hn = A.* Beta.^nn .* exp(-(w-(nn+1).*F).^2./(2*B));
                out = out + Hn;
            end
            
            out = out.*(1/A).*d2A_dAbar2;
        end
        
        function out=gaussian_harmonic_d2H_dAbar_dBbar(obj, w, X, ~)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            Beta = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dBeta_dBetabar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2Beta_dBetabar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Get number of harmonics
            N = obj.modelOptions;
            
            H0 = A.*exp(-(w-F).^2./(2*B));
            out = H0.*( ((w-F).^2)./(2*B^2) .* (1/A) );
            for nn = 1:N
                Hn = A.* Beta.^nn .* exp(-(w-(nn+1).*F).^2./(2*B));
                out = out + Hn.*( ((w-(nn+1)*F).^2)./(2*B.^2) .* (1/A) );
            end
            
            out = out.*dA_dAbar.*dB_dBbar;
        end
        
        function out=gaussian_harmonic_d2H_dAbar_dBetabar(obj, w, X, ~)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            Beta = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dBeta_dBetabar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2Beta_dBetabar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Get number of harmonics
            N = obj.modelOptions;
            
            H0 = A.*exp(-(w-F).^2./(2*B));
            fact = 0;
            out = H0.*fact;
            for nn = 1:N
                Hn = A.* Beta.^nn .* exp(-(w-(nn+1).*F).^2./(2*B));
                fact = (nn/Beta) .* (1/A);
                out = out + Hn.*fact;
            end
            
            out = out.*dBeta_dBetabar.*dA_dAbar;
        end
        
        function out=gaussian_harmonic_d2H_dBbar2(obj, w, X, ~)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            Beta = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dBeta_dBetabar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2Beta_dBetabar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Get number of harmonics
            N = obj.modelOptions;
            
            H0 = A.*exp(-(w-F).^2./(2*B));
            Z0 = ((w-F).^2)./(2*B.^2);
            dZ0_dB = -((w-F).^2)./(B.^3);
            first_fact = H0.*(Z0.^2 + dZ0_dB);
            dH_dB = H0.*Z0;
            for nn = 1:N
                Hn = A.* Beta.^nn .* exp(-(w-(nn+1).*F).^2./(2*B));
                Zn = ((w-(nn+1)*F).^2)./(2*B.^2);
                dZn_dB = -((w-(nn+1)*F).^2)./(B.^3);
                first_fact = first_fact + Hn.*(Zn.^2 + dZn_dB);
                dH_dB = dH_dB + Hn.*Zn;
            end
            
            out = first_fact.*(dB_dBbar.^2) + dH_dB.*d2B_dBbar2;
        end
        
        function out=gaussian_harmonic_d2H_dBbar_dBetabar(obj, w, X, ~)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            Beta = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dBeta_dBetabar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2Beta_dBetabar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Get number of harmonics
            N = obj.modelOptions;
            
            H0 = A.*exp(-(w-F).^2./(2*B));
            fact = 0;
            out = H0.*fact;
            for nn = 1:N
                Hn = A.* Beta.^nn .* exp(-(w-(nn+1).*F).^2./(2*B));
                fact = ((w-(nn+1)*F).^2)/(2*B.^2) .* (nn/Beta);
                out = out + Hn.*fact;
            end
            
            out = out.*dB_dBbar.*dBeta_dBetabar;
        end
        
        function out=gaussian_harmonic_d2H_dBetabar2(obj, w, X, ~)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            Beta = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dBeta_dBetabar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2Beta_dBetabar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Get number of harmonics
            N = obj.modelOptions;
            
            first_fact = 0;
            dH_dBeta = 0;
            for nn = 1:N
                Hn = A.* Beta.^nn .* exp(-(w-(nn+1).*F).^2./(2*B));
                Zn = nn/Beta;
                dZn_dBeta = -nn/Beta.^2;
                first_fact = first_fact + Hn.*(Zn.^2 + dZn_dBeta);
                dH_dBeta = dH_dBeta + Hn.*Zn;
            end
            
            out = first_fact.*dBeta_dBetabar.^2 + dH_dBeta.*d2Beta_dBetabar2;
        end
        
        
        %% BOX FUNCTION MODEL
        %-----------------------------------
        % EXPONENTIAL BOX FUNCTION
        %-----------------------------------
        function H_peak = boxfun(obj, w, X)
            if size(X,1) ~= 3
                error('Box model requires 3 parameters: F, A, and B');
            end
            
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            
            %Get number of harmonics
            p = obj.modelOptions;
            
            %Check if even
            if mod(p,2)
                error('Box function exponent must be even');
            end
            
            H_peak=A.*exp(-(w-F).^p./(2*B));
        end
        
        %--------------------------------------
        % EXPONENTIAL BOX FUNCTION DERIVATIVES
        %--------------------------------------
        function out=boxfun_dH_dFbar(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            
            %Get number of harmonics
            p = obj.modelOptions;
            
            out=H_peak.*p.*(w-F).^(p-1)/(2.*B).*dF_dFbar;
        end
        
        function out=boxfun_dH_dAbar(obj, ~, X, H_peak)
            %Get param values through link functions
            A = obj.modelParams(2).linkFun(X(2,:));
            
            %Get derivatives of link functions
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            
            out=H_peak./A.*dA_dAbar;
        end
        
        function out=boxfun_dH_dBbar(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get number of harmonics
            p = obj.modelOptions;
            
            out=H_peak.*(w-F).^p./(2*B.^2).*dB_dBbar;
        end
        
        
        %-----------------------------------
        %  EXPONENTIAL BOX MODEL 2nd DERIVAIVES
        %-----------------------------------
        function out = box_d2H_dFbar2(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            
            %Get number of harmonics
            p = obj.modelOptions;
            
            Z = p.*(w-F).^(p-1)./(2.*B);
            dZ_dF = - (p.*(p-1).*(w-F).^(p-2))./(2.*B);
            first_part =  (Z).^2  + dZ_dF;
            out = H_peak.*(first_part.*(dF_dFbar.^2) + Z.*d2F_dFbar2);
        end
        
        function out = box_d2H_dFbar_dAbar(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            
            %Get number of harmonics
            p = obj.modelOptions;
            
            Z = p*(w-F).^(p-1)./(2*B);
            out = H_peak .* Z ./ A .*dF_dFbar .* dA_dAbar;
        end
        
        function out = box_d2H_dFbar_dBbar(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            
            %Get number of harmonics
            p = obj.modelOptions;
            
            out = H_peak .*((p*(w-F).^(2*p-1))./(4*B.^3) - (p*(w-F).^(p-1))./(2*B.^2)) .*dF_dFbar .*dB_dBbar;
        end
        
        function out = box_d2H_dAbar2(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            
            %Get number of harmonics
            p = obj.modelOptions;
            
            out = H_peak./A.*d2A_dAbar2;
        end
        
        function out = box_d2H_dAbar_dBbar(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            
            %Get number of harmonics
            p = obj.modelOptions;
            
            out = H_peak.*( ((w-F).^p)./(2*B.^2) .* (1/A) ) .* dA_dAbar .* dB_dBbar;
        end
        
        function out = box_d2H_dBbar2(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            
            %Get number of harmonics
            p = obj.modelOptions;
            
            Z = ((w-F).^p)./(2*B.^2);
            dZ_dB = -((w-F).^p)./(B.^3);
            
            out = H_peak.*((Z.^2 + dZ_dB).*(dB_dBbar.^2) + Z.*d2B_dBbar2);
        end
        
        
        %% SKEWED GAMMA MODEL
        
        %-----------------------------------
        %  SKEWED GAMMA FUNCITON
        %-----------------------------------
        function H_peak = skewgamma(obj, w, X)
            if size(X,1) ~= 4
                error('Skewed gamma model requires 4 parameters: F, A, B and S');
            end
            
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            S = obj.modelParams(4).linkFun(X(4,:));
            
            %Compute the gamma parameters
            Alpha = 4./S.^2;
            Beta = 2./(S.*sqrt(B));
            O = F + S.*sqrt(B)./2 - 2.*sqrt(B)./S;
            
            %Compute the normalization factor
            N = (Alpha-1).*(log(Alpha-1)-log(Beta)-1);
            
            %Compute the function
            H_peak = A.*exp((Alpha-1).*log(w-O)-Beta.*(w-O)-N);
            %Set all pre-shift values to zero
            H_peak(w<=O) = 0;
        end
        
        %-----------------------------------
        % SKEWED GAMMA FUNCITON DERIVATIVES
        %-----------------------------------
        %-----------------------------------
        %  dH/_dFbar
        %-----------------------------------
        function out=skewgamma_dH_dFbar(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            S = obj.modelParams(4).linkFun(X(4,:));
            
            %Compute the gamma parameters
            Alpha = 4./S.^2;
            Beta = 2./(S.*sqrt(B));
            O = F + S.*sqrt(B)./2 - 2.*sqrt(B)./S;
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            
            dH_dO = H_peak.*(-(Alpha-1)./(w-O)+Beta);
            out = dH_dO*dF_dFbar;
            
            if abs(out)>1e99
                warning('dH/_dFbar blown up');
            end
        end
        
        %-----------------------------------
        %  dH/_dAbar
        %-----------------------------------
        function out=skewgamma_dH_dAbar(obj, ~, X, H_peak)
            %Get param values through link functions
            A = obj.modelParams(2).linkFun(X(2,:));
            
            %Get derivatives of link functions
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            
            out=H_peak./A * dA_dAbar;
            
            if abs(out)>1e99
                warning('dH/_dAbar blown up');
            end
        end
        
        %-----------------------------------
        %  dH/_dBbar
        %-----------------------------------
        function out=skewgamma_dH_dBbar(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            S = obj.modelParams(4).linkFun(X(4,:));
            
            %Compute the gamma parameters
            Alpha = 4./S.^2;
            Beta = 2./(S.*sqrt(B));
            O = F + S.*sqrt(B)./2 - 2.*sqrt(B)./S;
            
            %Get derivatives of link functions
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            dH_dBeta = H_peak.*(-(w-O) +(Alpha-1)./Beta);
            dBeta_dB = -Beta./(2*B);
            dH_dO = H_peak.*(-(Alpha-1)./(w-O)+Beta);
            dO_dB = -1./(S.*sqrt(B))+S./(4.*sqrt(B));
            
            out = (dH_dBeta.*dBeta_dB + dH_dO.*dO_dB).*dB_dBbar;
            
            if abs(out)>1e99
                warning('dH/_dBbar blown up');
            end
        end
        
        
        %-----------------------------------
        %  dH/_dSbar
        %-----------------------------------
        function out=skewgamma_dH_dSbar(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            S = obj.modelParams(4).linkFun(X(4,:));
            
            %Compute the gamma parameters
            Alpha = 4./S.^2;
            Beta = 2./(S.*sqrt(B));
            O = F + S.*sqrt(B)./2 - 2.*sqrt(B)./S;
            
            %Get derivatives of link functions
            dS_dSbar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            dH_dAlpha = H_peak.*(log(w-O) - log((Alpha-1)/Beta));
            dAlpha_dS = -8./S.^3;
            dH_dBeta = H_peak.*(-(w-O) +(Alpha-1)./Beta);
            dBeta_dS = -Beta/S;
            dH_dO = H_peak.*(-(Alpha-1)./(w-O)+Beta);
            dO_dS = (2.*sqrt(B))./(S.^2)+sqrt(B)./2;
            
            out=(dH_dAlpha.*dAlpha_dS + dH_dBeta.*dBeta_dS + dH_dO.*dO_dS).*dS_dSbar;
            
            if abs(out)>1e99
                warning('dH/_dSbar blown up');
            end
        end
        
        %-----------------------------------
        % SKEWED GAMMA FUNCITON 2nd DERIVATIVES
        %-----------------------------------
        function out = skewgamma_d2H_dFbar2(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            S = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dS_dSbar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2S_dSbar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Compute the gamma parameters
            Alpha = 4./S.^2;
            Beta = 2./(S.*sqrt(B));
            O = F + S.*sqrt(B)./2 - 2.*sqrt(B)./S;
            
            Z = Beta-(Alpha-1)./(w-O);
            dZ_dF = -(Alpha-1)./((w-O).^2);
            
            out = H_peak.*( (Z.^2 + dZ_dF).*(dF_dFbar.^2) + Z.*d2F_dFbar2 );
        end
        
        function out = skewgamma_d2H_dFbar_dAbar(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            S = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dS_dSbar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2S_dSbar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Compute the gamma parameters
            Alpha = 4./S.^2;
            Beta = 2./(S.*sqrt(B));
            O = F + S.*sqrt(B)./2 - 2.*sqrt(B)./S;
            
            out = H_peak.*(Beta-(Alpha-1)./(w-O)).*(1/A).*dF_dFbar.*dA_dAbar;
            
        end
        
        function out = skewgamma_d2H_dFbar_dBbar(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            S = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dS_dSbar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2S_dSbar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Compute the gamma parameters
            Alpha = 4./S.^2;
            Beta = 2./(S.*sqrt(B));
            O = F + S.*sqrt(B)./2 - 2.*sqrt(B)./S;
            
            Theta = (Beta-(Alpha-1)./(w-O));
            Psi = ((w-O).*Beta.^3)./(2.*Alpha) -((Alpha-1).*Beta.^2)./(Alpha) +(Beta.*(Alpha-1).^2)./(2*Alpha.*(w-O));
            
            out = H_peak.*( Psi.*Theta +(Beta.*(Alpha-1).^2)./(2*Alpha.*(w-O).^2) - (Beta.^3)./(2*Alpha) ) .*dF_dFbar.*dB_dBbar;
        end
        
        function out = skewgamma_d2H_dFbar_dSbar(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            S = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dS_dSbar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2S_dSbar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Compute the gamma parameters
            Alpha = 4./S.^2;
            Beta = 2./(S.*sqrt(B));
            O = F + S.*sqrt(B)./2 - 2.*sqrt(B)./S;
            
            Theta = (Beta-(Alpha-1)./(w-O));
            Delta = -Alpha.^(3/2).*log((Beta.*(w-O))./(Alpha-1))+Alpha.^(1/2)+Beta.*(Alpha.^(1/2)).*(w-O)/2-(Alpha.^(5/2)-Alpha.^(1/2))./(2*Beta.*(w-O));
            
            out = H_peak.*(Delta.*Theta +(Alpha.^(3/2))./(w-O) -(Beta.*Alpha.^(1/2))/2 -(Alpha.^(5/2)-Alpha.^(1/2))./(2*Beta.*(w-O).^2)).*dF_dFbar.*dS_dSbar;
        end
        
        function out = skewgamma_d2H_dAbar2(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            S = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dS_dSbar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2S_dSbar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Compute the gamma parameters
            Alpha = 4./S.^2;
            Beta = 2./(S.*sqrt(B));
            O = F + S.*sqrt(B)./2 - 2.*sqrt(B)./S;
            
            out=H_peak./A * d2A_dAbar2;
        end
        
        function out = skewgamma_d2H_dAbar_dBbar(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            S = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dS_dSbar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2S_dSbar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Compute the gamma parameters
            Alpha = 4./S.^2;
            Beta = 2./(S.*sqrt(B));
            O = F + S.*sqrt(B)./2 - 2.*sqrt(B)./S;
            
            Psi = ((w-O).*Beta.^3)./(2.*Alpha)-((Alpha-1).*Beta.^2)./(Alpha) + (Beta.*(Alpha-1).^2)./(2*Alpha.*(w-O));
            
            out = H_peak.*( Psi./A ) .*dA_dAbar.*dB_dBbar;
        end
        
        function out = skewgamma_d2H_dAbar_dSbar(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            S = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dS_dSbar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2S_dSbar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Compute the gamma parameters
            Alpha = 4./S.^2;
            Beta = 2./(S.*sqrt(B));
            O = F + S.*sqrt(B)./2 - 2.*sqrt(B)./S;
            
            Delta = -Alpha.^(3/2).*log((Beta.*(w-O))./(Alpha-1))+Alpha.^(1/2)+Beta.*(Alpha.^(1/2)).*(w-O)/2-(Alpha.^(5/2)-Alpha.^(1/2))./(2*Beta.*(w-O));
            
            out = H_peak.*(Delta./A).*dA_dAbar.*dS_dSbar;
        end
        
        function out = skewgamma_d2H_dBbar2(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            S = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dS_dSbar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2S_dSbar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Compute the gamma parameters
            Alpha = 4./S.^2;
            Beta = 2./(S.*sqrt(B));
            O = F + S.*sqrt(B)./2 - 2.*sqrt(B)./S;
            
            Psi = ((w-O).*Beta.^3)./(2.*Alpha) -((Alpha-1).*Beta.^2)./(Alpha) +(Beta.*(Alpha-1).^2)./(2*Alpha.*(w-O));
            
            dPsi_dB= -((Beta.^2 .* (-1 + Alpha + Beta .* (O - w)).^2 .* (-1 + Alpha + 3 .* Beta.* (-O + w)))./(4 .* Alpha.^2 .* (O - w).^2));
            first_part = Psi.^2 + dPsi_dB;
            
            out = H_peak.*( first_part.*(dB_dBbar.^2) + Psi.*d2B_dBbar2 );
        end
        
        
        function out = skewgamma_d2H_dBbar_dSbar(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            S = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dS_dSbar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2S_dSbar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Compute the gamma parameters
            Alpha = 4./S.^2;
            Beta = 2./(S.*sqrt(B));
            O = F + S.*sqrt(B)./2 - 2.*sqrt(B)./S;
            
            Psi = ((w-O).*Beta.^3)./(2.*Alpha)-((Alpha-1).*Beta.^2)./(Alpha) + (Beta.*(Alpha-1).^2)./(2*Alpha.*(w-O));
            Delta = -Alpha.^(3/2).*log((Beta.*(w-O))./(Alpha-1))+Alpha.^(1/2)+Beta.*(Alpha.^(1/2)).*(w-O)/2-(Alpha.^(5/2)-Alpha.^(1/2))./(2*Beta.*(w-O));
            
            fact = Psi.*Delta -(Beta.^3.*(w-O))/(4.*sqrt(Alpha))-(Beta.*(3*Alpha+1).*(Alpha-1))./(4*sqrt(Alpha).*(w-O)) ...
                +(Beta.^2.*(3*Alpha-1))/(4*sqrt(Alpha)) +((Alpha-1).^2.*(Alpha+1))./(4*sqrt(Alpha).*(w-O).^2);
            
            out = H_peak.*fact.*dB_dBbar.*dS_dSbar;
        end
        
        function out = skewgamma_d2H_dSbar2(obj, w, X, H_peak)
            %Get param values through link functions
            F = obj.modelParams(1).linkFun(X(1,:));
            A = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            S = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dF_dFbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dA_dAbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dB_dBbar = obj.modelParams(3).linkDerivFun(X(3,:));
            dS_dSbar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Get 2nd derivatives of link functions
            d2F_dFbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2A_dAbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2B_dBbar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            d2S_dSbar2 = obj.modelParams(4).linkDeriv2Fun(X(4,:));
            
            %Compute the gamma parameters
            Alpha = 4./S.^2;
            Beta = 2./(S.*sqrt(B));
            O = F + S.*sqrt(B)./2 - 2.*sqrt(B)./S;
            
            Delta = -Alpha.^(3/2).*log((Beta.*(w-O))./(Alpha-1))+Alpha.^(1/2)+Beta.*(Alpha.^(1/2)).*(w-O)/2-(Alpha.^(5/2)-Alpha.^(1/2))./(2*Beta.*(w-O));
            
            first_part = Delta.^2 +(3*Alpha.^2)/2.*log((Beta.*(w-O))./(Alpha-1)) -(Alpha.*(3*Alpha.^2+4*Alpha-3))./(4*(Alpha-1)) ...
                -(Beta.*Alpha.*(w-O))/2 +(Alpha.^2.*(3*Alpha+1))./(2*Beta.*(w-O)) -(Alpha.*(Alpha.^2-1).*(Alpha+1))./(4*Beta.^2.*(w-O).^2);
            
            out = H_peak.*(first_part.*(dS_dSbar.^2) + Delta.*d2S_dSbar2);
        end
        
        %% LINEAR MODEL
        
        %-----------------------------------
        %  LINEAR
        %-----------------------------------
        function H_peak = linear(obj, w, X, ~)
            if size(X,1) ~= 2
                error('Linear model requires 2 parameters: M and B');
            end
            
            %Get param values through link functions
            M = obj.modelParams(1).linkFun(X(1,:));
            B = obj.modelParams(2).linkFun(X(2,:));
            
            %Compute the function
            H_peak = M.*w +B;
        end
        
        %-----------------------------------
        %  LINEAR DERIVATIVES
        %-----------------------------------
        function out = linear_dH_dMbar(obj, w, X, ~)
            %Get derivatives of link functions
            dM_dMbar = obj.modelParams(1).linkDerivFun(X(1,:));
            
            out = w.*dM_dMbar;
        end
        
        function out = linear_dH_dBbar(~, w, ~, ~)
            out = ones(size(w));
        end
        
        %-----------------------------------
        %  LINEAR 2nd DERIVATIVES
        %-----------------------------------
        function out = linear_dH_dMbar2(~, w, ~, ~)
            out = zeros(size(w));
        end
        
        function out = linear_dH_dMbar_dBbar(~, w, ~, ~)
            out = zeros(size(w));
        end
        
        function out = linear_dH_dBbar2(~, w, ~, ~)
            out = zeros(size(w));
        end
        
        %% EXPONENTIAL DECAY MODEL
        
        %-----------------------------------
        %  EXPONENIAL DECAY
        %-----------------------------------
        function H_peak = expdecay(obj, w, X)
            if size(X,1) ~= 3
                error('Exponential decay model requires 3 parameters: BA, BR, and BO');
            end
            
            %Get param values through link functions
            BA = obj.modelParams(1).linkFun(X(1,:));
            BR = obj.modelParams(2).linkFun(X(2,:));
            BO = obj.modelParams(3).linkFun(X(3,:));
            
            %Compute the function
            H_peak = BA.*(1-BR).^w +BO;
        end
        
        %-----------------------------------
        %  EXPONENIAL DECAY DERIVATIVES
        %-----------------------------------
        function out = expdecay_dHdBA_bar(obj, w, X, ~)
            %Get param values through link functions
            BR = obj.modelParams(2).linkFun(X(2,:));
            
            %Get derivatives of link functions
            dBA_dBAbar = obj.modelParams(1).linkDerivFun(X(1,:));
            
            out=(1-BR).^w .* dBA_dBAbar;
        end
        
        function out = expdecay_dHdBR_bar(obj, w, X, ~)
            %Get param values through link functions
            BA = obj.modelParams(1).linkFun(X(1,:));
            BR = obj.modelParams(2).linkFun(X(2,:));
            
            %Get derivatives of link functions
            dBR_dBRbar = obj.modelParams(2).linkDerivFun(X(2,:));
            
            out = -BA.*w.*(1-BR).^(w-1).*dBR_dBRbar;
        end
        
        function out = expdecay_dHdBO_bar(obj, w, X, ~)
            %Get derivatives of link functions
            dBO_dBObar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            out = ones(size(w)).*dBO_dBObar;
        end
        
        %-----------------------------------
        %  EXPONENIAL DECAY 2nd DERIVATIVES
        %-----------------------------------
        function out = expdecay_d2H_dBAbar2(obj, w, X, ~)
            %Get param values through link functions
            BA = obj.modelParams(1).linkFun(X(1,:));
            BR = obj.modelParams(2).linkFun(X(2,:));
            BO = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dBA_dBAbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dBR_dBRbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dBO_dBObar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get 2nd derivatives of link functions
            d2BA_dBAbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2BR_dBRbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2BO_dBObar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            
            out = (1-BR).^w .* d2BA_dBAbar2;
        end
        
        function out = expdecay_d2H_dBAbar_dBRbar(obj, w, X, ~)
            %Get param values through link functions
            BA = obj.modelParams(1).linkFun(X(1,:));
            BR = obj.modelParams(2).linkFun(X(2,:));
            BO = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dBA_dBAbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dBR_dBRbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dBO_dBObar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get 2nd derivatives of link functions
            d2BA_dBAbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2BR_dBRbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2BO_dBObar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            
            out = -w.*(1-BR).^(w-1) .* dBA_dBAbar .*dBR_dBRbar;
        end
        
        function out = expdecay_d2H_dBAbar_dBObar(obj, w, X, ~)
            out = zeros(size(w));
        end
        
        function out = expdecay_d2H_dBRbar2(obj, w, X, ~)
            %Get param values through link functions
            BA = obj.modelParams(1).linkFun(X(1,:));
            BR = obj.modelParams(2).linkFun(X(2,:));
            BO = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dBA_dBAbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dBR_dBRbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dBO_dBObar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get 2nd derivatives of link functions
            d2BA_dBAbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2BR_dBRbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2BO_dBObar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            
            out = BA.*w.*(w-1).*(1-BR).^(w-2).*(-dBR_dBRbar).^2  -BA.*w.*(1-BR).^(w-1).*d2BR_dBRbar2;
        end
        
        function out = expdecay_d2H_dBRbar_dBObar(~, w, ~, ~)
            out = zeros(size(w));
        end
        
        function out = expdecay_d2H_dBObar2(obj, w, X, ~)
            %Get param values through link functions
            BA = obj.modelParams(1).linkFun(X(1,:));
            BR = obj.modelParams(2).linkFun(X(2,:));
            BO = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dBA_dBAbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dBR_dBRbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dBO_dBObar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get 2nd derivatives of link functions
            d2BA_dBAbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2BR_dBRbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2BO_dBObar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            
            out = ones(size(w)).*d2BO_dBObar2;
        end
        
        %% -----------------------------------
        %  DOUBLE EXPONENIAL DEAY
        %-----------------------------------
        function H_peak = doubledecay(obj, w, X)
            if size(X,1) ~= 4
                error('Double exponential decay model requires 5 parameters: BA1, BR1, BA2, BR2, and BO');
            end
            
            %Get param values through link functions
            BA1 = obj.modelParams(1).linkFun(X(1,:));
            BR1 = obj.modelParams(2).linkFun(X(2,:));
            BA2 = obj.modelParams(3).linkFun(X(3,:));
            BR2 = obj.modelParams(4).linkFun(X(4,:));
            BO = obj.modelParams(5).linkFun(X(5,:));
            
            %Compute the function
            H_peak = BA1.*(1-BR1).^w + BA2.*(1-BR2).^w +BO;
        end
        
        %-----------------------------------
        %  DOUBLE EXPONENIAL DEAY DERIVATIVES
        %-----------------------------------
        function out = doubledecay_dHdBA1_bar(obj, w, X, ~)
            %Get param values through link functions
            BR1 = obj.modelParams(2).linkFun(X(2,:));
            
            %Get derivatives of link functions
            dBA1_dBA1bar = obj.modelParams(1).linkDerivFun(X(1,:));
            
            out = (1-BR1).^w .* dBA1_dBA1bar;
        end
        
        function out=doubledecay_dHdBR1_bar(obj, w, X, ~)
            %Get param values through link functions
            BA1 = obj.modelParams(1).linkFun(X(1,:));
            BR1 = obj.modelParams(2).linkFun(X(2,:));
            
            %Get derivatives of link functions
            dBR1_dBR1bar = obj.modelParams(2).linkDerivFun(X(2,:));
            
            out = -BA1.*w.*(1-BR1).^(w-1).*dBR1_dBR1bar;
        end
        
        function out=doubledecay_dHdBA2_bar(obj, w, X, ~)
            %Get param values through link functions
            BR2 = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dBA2_dBA2bar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            out = (1-BR2).^w .* dBA2_dBA2bar;
        end
        
        function out=doubledecay_dHdBR2_bar(obj, w, X, ~)
            %Get param values through link functions
            BA2 = obj.modelParams(3).linkFun(X(3,:));
            BR2 = obj.modelParams(4).linkFun(X(4,:));
            
            %Get derivatives of link functions
            dBR2_dBR2bar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            out = -BA2.*w.*(1-BR2).^(w-1).*dBR2_dBR2bar;
        end
        
        function out=doubledecay_dHdBO_bar(obj, w, X, ~)
            %Get derivatives of link functions
            dBO_dBObar = obj.modelParams(5).linkDerivFun(X(5,:));
            
            out = ones(size(w))*dBO_dBObar;
        end
        
        %-----------------------------------
        %  SMOOTHED DOUBLE EXPONENIAL
        %-----------------------------------
        function H_peak = smoothdoubledecay(obj, w, X)
            if size(X,1) ~= 4
                error('Gamma double exponential decay model requires 4 parameters: L1mL2, L2, K, and O');
            end
            
            %Get param values through link functions
            L1mL2 = obj.modelParams(1).linkFun(X(1,:));
            L2 = obj.modelParams(2).linkFun(X(2,:));
            K = obj.modelParams(3).linkFun(X(3,:));
            O = obj.modelParams(4).linkFun(X(4,:));
            
            %Set alpha and theta as fixed values
            Alpha = obj.modelOptions(1);
            Theta = obj.modelOptions(2);
            
            Beta=1./Theta;
            
            combined_rate = Beta + L1mL2;
            
            %Computes the weight function
            W1 = exp(-(L1mL2+L2).*w) .* gammainc(Beta.*w, Alpha ,'upper');
            if combined_rate==0
                W2 = exp(-L2.*w) .* w.^Alpha .* (Beta .^ Alpha) ./ gamma(Alpha+1);
            elseif combined_rate>0
                W2 = (Beta ./ combined_rate) .^ Alpha .* exp(-L2.*w) .* (1 - gammainc(combined_rate.*w, Alpha ,'upper'));
            else
                W2 = -(Beta ./ combined_rate) .^ Alpha .* exp(-L2.*w) .* (1 - gammainc(combined_rate.*w, Alpha ,'upper'));
            end
            
            
            %Full function
            H_peak = K.*(W1+W2)+O;
        end
        
        %-----------------------------------------------
        %  SMOOTHED DOUBLE EXPONENIAL DERIVAIVES
        %-----------------------------------------------
        function out = smoothdoubledecay_dHdL1mL2t_bar(obj, w, X, ~)
            %Get param values through link functions
            L1mL2 = obj.modelParams(1).linkFun(X(1,:));
            L2 = obj.modelParams(2).linkFun(X(2,:));
            K = obj.modelParams(3).linkFun(X(3,:));
            
            %Set alpha and theta as fixed values
            Alpha = obj.modelOptions(1);
            Theta = obj.modelOptions(2);
            
            %Get derivatives of link functions
            dL1m2_dL1m2bar = obj.modelParams(1).linkDerivFun(X(1,:));
            
            Beta = 1/Theta;
            
            combined_rate = Beta + L1mL2;
            
            %Computes the weight function
            W1 = exp(-(L1mL2+L2)*w) .* gammainc(Beta*w, Alpha ,'upper');
            if combined_rate==0
                W2 = exp(-L2*w) .* w^Alpha * (Beta ^ Alpha) / gamma(Alpha+1);
            elseif combined_rate>0
                W2 = (Beta / combined_rate) ^ Alpha * exp(-L2*w) .* (1 - gammainc(combined_rate*w, Alpha ,'upper'));
            else
                W2 = -(Beta / combined_rate) ^ Alpha * exp(-L2*w) .* (1 - gammainc(combined_rate*w, Alpha ,'upper'));
            end
            
            %Full function
            term1 = -w.*W1;
            if combined_rate==0
                term2 = 0;
                term3 = 0;
            else
                term2 = -Alpha*W2/combined_rate;
                term3 = exp(-(Beta + L1mL2)*w) .* w.^Alpha * Beta .^ Alpha / (gamma(Alpha) * combined_rate);
            end
            out = K *(term1 +term2 + term3) * dL1m2_dL1m2bar;
        end
        
        function out = smoothdoubledecay_dHdL2_bar(obj, w, X, ~)
            %Get param values through link functions
            L1mL2 = obj.modelParams(1).linkFun(X(1,:));
            L2 = obj.modelParams(2).linkFun(X(2,:));
            K = obj.modelParams(3).linkFun(X(3,:));
            
            %Set alpha and theta as fixed values
            Alpha = obj.modelOptions(1);
            Theta = obj.modelOptions(2);
            
            %Get derivatives of link functions
            dL2_dL2bar = obj.modelParams(2).linkDerivFun(X(2,:));
            
            Beta = 1/Theta;
            
            combined_rate = Beta + L1mL2;
            
            %Computes the weight function
            W1 = exp(-(L1mL2+L2)*w) .* gammainc(Beta*w, Alpha ,'upper');
            if combined_rate==0
                W2 = exp(-L2*w) .* w.^Alpha * (Beta .^ Alpha) / gamma(Alpha+1);
            elseif combined_rate>0
                W2 = (Beta / combined_rate) ^ Alpha * exp(-L2*w) .* (1 - gammainc(combined_rate*w, Alpha, 'upper'));
            else
                W2 = -(Beta / combined_rate) ^ Alpha * exp(-L2*w) .* (1 - gammainc(combined_rate*w, Alpha, 'upper'));
            end
            
            %Full function
            term1 = -w.*W1;
            term2 = -w.*W2;
            out = K*(term1 +term2) * dL2_dL2bar;
        end
        
        function out = smoothdoubledecay_dHdK_bar(obj, w, X, ~)
            %Get param values through link functions
            L1mL2 = obj.modelParams(1).linkFun(X(1,:));
            L2 = obj.modelParams(2).linkFun(X(2,:));
            
            %Set alpha and theta as fixed values
            Alpha = obj.modelOptions(1);
            Theta = obj.modelOptions(2);
            
            Beta = 1/Theta;
            
            %Get derivatives of link functions
            dK_dKbar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            combined_rate = Beta + L1mL2;
            
            %Computes the weight function
            W1 = exp(-(L1mL2+L2).*w) .* gammainc(Beta.*w, Alpha ,'upper');
            
            if combined_rate==0
                W2 = exp(-L2.*w) .* w.^Alpha .* (Beta .^ Alpha) ./ gamma(Alpha+1);
            elseif combined_rate>0
                W2 = (Beta ./ combined_rate) .^ Alpha .* exp(-L2.*w) .* (1 - gammainc(combined_rate.*w, Alpha ,'upper'));
            else
                W2 = -(Beta ./ combined_rate) .^ Alpha .* exp(-L2.*w) .* (1 - gammainc(combined_rate.*w, Alpha ,'upper'));
            end
            
            
            %Full function
            out = (W1+W2)*dK_dKbar;
        end
        
        function out = smoothdoubledecay_dHdO_bar(obj, w, X, ~)
            %Get derivatives of link functions
            dO_dObar = obj.modelParams(4).linkDerivFun(X(4,:));
            
            %Full function
            out = ones(size(w))*dO_dObar;
        end
        
        %% 1/F MODEL
        %-----------------------------------
        %  1/F DECAY
        %-----------------------------------
        function H_peak = oof(obj, w, X)
            if size(X,1) ~= 3
                error('One-over-F model requires 3 parameters: K, R, and B');
            end
            
            %Get param values through link functions
            K = obj.modelParams(1).linkFun(X(1,:));
            R = obj.modelParams(2).linkFun(X(2,:));
            B = obj.modelParams(3).linkFun(X(3,:));
            
            %Compute the function
            H_peak = B-log(K+w.^R);
        end
        
        %-----------------------------------
        %  1/F DERIVATIVES
        %-----------------------------------
        function out = oof_dHdK_bar(obj, w, X, ~)
            %Get param values through link functions
            K = obj.modelParams(1).linkFun(X(1,:));
            R = obj.modelParams(2).linkFun(X(2,:));
            
            %Get derivatives of link functions
            dK_dKbar = obj.modelParams(1).linkDerivFun(X(1,:));
            
            out = -1./(K+w.^R) .* dK_dKbar;
        end
        
        function out = oof_dHdR_bar(obj, w, X, ~)
            %Get param values through link functions
            K = obj.modelParams(1).linkFun(X(1,:));
            R = obj.modelParams(2).linkFun(X(2,:));
            
            %Get derivatives of link functions
            dR_dRbar = obj.modelParams(2).linkDerivFun(X(2,:));
            
            out = -1./(K+w.^R) .* (w.^R) .* log(w).* dR_dRbar;
        end
        
        function out = oof_dHdB_bar(obj, w, X, ~)
            %Get derivatives of link functions
            dO_dObar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            out = ones(size(w))*dO_dObar;
        end
        
        %-----------------------------------
        %  1/F 2nd DERIVATIVES
        %-----------------------------------
        function out = oof_d2H_dKbar2(obj, w, X, ~)
            %Get param values through link functions
            BA = obj.modelParams(1).linkFun(X(1,:));
            BR = obj.modelParams(2).linkFun(X(2,:));
            BO = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dBA_dBAbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dBR_dBRbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dBO_dBObar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get 2nd derivatives of link functions
            d2BA_dBAbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2BR_dBRbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2BO_dBObar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            
            out = (1-BR).^w .* d2BA_dBAbar2;
        end
        
        function out = oof_d2H_dKbar_dRbar(obj, w, X, ~)
            %Get param values through link functions
            BA = obj.modelParams(1).linkFun(X(1,:));
            BR = obj.modelParams(2).linkFun(X(2,:));
            BO = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dBA_dBAbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dBR_dBRbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dBO_dBObar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get 2nd derivatives of link functions
            d2BA_dBAbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2BR_dBRbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2BO_dBObar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            
            out = -w.*(1-BR).^(w-1) .* dBA_dBAbar .*dBR_dBRbar;
        end
        
        function out = oof_d2H_dKbar_dBbar(~, w, ~, ~)
            out = zeros(size(w));
        end
        
        function out = oof_d2H_dRbar2(obj, w, X, ~)
            %Get param values through link functions
            BA = obj.modelParams(1).linkFun(X(1,:));
            BR = obj.modelParams(2).linkFun(X(2,:));
            BO = obj.modelParams(3).linkFun(X(3,:));
            
            %Get derivatives of link functions
            dBA_dBAbar = obj.modelParams(1).linkDerivFun(X(1,:));
            dBR_dBRbar = obj.modelParams(2).linkDerivFun(X(2,:));
            dBO_dBObar = obj.modelParams(3).linkDerivFun(X(3,:));
            
            %Get 2nd derivatives of link functions
            d2BA_dBAbar2 = obj.modelParams(1).linkDeriv2Fun(X(1,:));
            d2BR_dBRbar2 = obj.modelParams(2).linkDeriv2Fun(X(2,:));
            d2BO_dBObar2 = obj.modelParams(3).linkDeriv2Fun(X(3,:));
            
            out = BA.*w.*(w-1).*(1-BR).^(w-2).*(-dBR_dBRbar).^2  -BA.*w.*(1-BR).^(w-1).*d2BR_dBRbar2;
        end
        
        function out = oof_d2H_dRbar_dBbar(~, w, ~, ~)
            out = zeros(size(w));
        end
        
        function out = oof_d2H_dBbar2(~, w, ~, ~)
            out = zeros(size(w));
        end
        
    end
end

