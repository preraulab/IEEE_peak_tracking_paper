classdef MultiPeakModel < handle
    %MultiPeakModel Parametric model comprised of the summation of multiple
    %PeakModel
    %
    %   obj = MultiPeakModel(peakModels)
    %
    %   Input:
    %       peakModels: array of class PeakModel
    %
    %   Copyright 2019 Prerau Laboratory
    %
    %   Last modified 5/17/2019 by Michael Prerau
    
    properties (SetAccess = private)
        numPeaks = 0; %Number of peaks in the model (cannot set)
        totalPeakParams = 0;  %Number of total parameters in full model across all peaks (cannot set)
    end
    
    properties (Access = public)
        peakModels = []; %Vector of PeakModels
    end
    
    methods (Access = public)
        %-----------------------------------
        %       CONSTRUCTOR
        %-----------------------------------
        function obj = MultiPeakModel(peakModels)
            %Constructor for the MultiPeakModel class
            %
            %   obj = MultiPeakModel(peakModels)
            %
            %   Input:
            %       peakModels: 1 x P vector of PeakModel, the new peaks to add
            
            if nargin == 0
                obj.peakModels = [];
                return;
            end
            
            for ii = 1:length(peakModels)
                obj.addPeakModel(peakModels(ii));
            end
        end
        
        %-----------------------------------
        %       ADD A PEAK TO THE MODEL
        %-----------------------------------
        function obj = addPeakModel(obj, p)
            %Adds a PeakModel to the MultiPeakModel class
            %
            %   obj = addPeakModel(newPeak)
            %
            %   Input:
            %       newPeak: PeakModel, the new peak to add
            
            
            obj.numPeaks = obj.numPeaks + 1;
            obj.totalPeakParams = obj.totalPeakParams + p.numPeakParams;
            
            if isempty(obj.peakModels)
                obj.peakModels = p;
            else
                obj.peakModels(obj.numPeaks) = p;
            end
        end
        
        %-----------------------------------
        % GET X INDICES RELATED TO A PEAK
        %-----------------------------------
        function idxs = getPeakIdxs(obj, pickPeaks)
            %Get the indices for the states in X associated with a given peak
            %
            %   idxs = getPeakIdxs(pickPeaks)
            %
            %   Input:
            %       pickPeaks: integer, the index of the peaks to add or
            %       logical selector
            %
            %   Output:
            %       idxs: 1 x totalPeakParams logical, index vector to X
            
            %Allocate the logical to false
            idxs = false(1, obj.totalPeakParams);
            
            if isempty(pickPeaks) || ~any(pickPeaks)
                return;
            end
            
            if islogical(pickPeaks) || isempty(setdiff(unique(pickPeaks),[0 1]))
                pickPeaks = find(pickPeaks);
                
                if max(pickPeaks) > obj.numPeaks || min(pickPeaks) <= 0
                    error('Invalid peak selected');
                end
            end
            
            for pp = 1:length(pickPeaks)
                %Find the start and end indices for the peak parameters
                if pickPeaks(pp) == 1
                    idx_start = 1;
                else
                    idx_start = sum([obj.peakModels(1:pickPeaks(pp)-1).numPeakParams]) + 1;
                end
                
                idx_end = sum([obj.peakModels(1:pickPeaks(pp)).numPeakParams]);
                
                %Create a logical for those indicies
                idxs(idx_start:idx_end) = true;
            end
        end
        
        %-----------------------------------
        % GET PARAMETER CLASSES
        %-----------------------------------
        function classes = getParamClasses(obj)
            classes = cell(obj.totalPeakParams,1);
            count = 0;
            for pp = 1:obj.numPeaks
                classes((count+1):(count+obj.peakModels(pp).numPeakParams)) = {obj.peakModels(pp).modelParams.linkFunClass};
                count = count + obj.peakModels(pp).numPeakParams;
            end
        end
        
        %-----------------------------------
        % GET PARAMETER BOUNDS
        %-----------------------------------
        function bounds = getParamBounds(obj)
            bounds = zeros(obj.totalPeakParams,2);
            count = 0;
            for pp = 1:obj.numPeaks
                for qq = 1:obj.peakModels(pp).numPeakParams
                    tmp = obj.peakModels(pp).modelParams(qq).linkFunParams;
                    if ~isempty(tmp)
                        bounds(count+1,:) = tmp;
                    end
                    count = count + 1;
                end
            end
        end
        
        
           %-----------------------------------
        % GET PARAMETER NAMES
        %-----------------------------------
        function paramNames = getParamNames(obj)
            paramNames = cat(2,obj.peakModels.paramNames);
        end
        
        
                %-----------------------------------
        % GET PARAMETER DESCRIPTIONS
        %-----------------------------------
        function paramDescriptions = getParamDescriptions(obj)
            paramDescriptions = cat(2,obj.peakModels.paramDescriptions);
        end
   
   
        
        %-----------------------------------
        % MODEL FUNCTION DERIVATIVE (H)
        %-----------------------------------
        function [H, components] = getH(obj, X, w, combo)
            %Compute the multipeak function and associated components
            %
            %   [H, components] = getH(X, w, combo)
            %
            %   Input:
            %       X: 1 x totalPeakParams double, state vector
            %       w: 1 x W double, frequency vector
            %       combo: 1 x numPeaks logical, selector for active peaks
            %
            %   Output:
            %       H: 1 x w double, multipeak function value
            %       components: numPeaks x w double, components from each
            %       of the individual peaks
            
            %Set default to all combos being true
            if nargin<4 || isempty(combo)
                combo = true(1,obj.totalPeakParams);
            end
            
            if nargin<3
                w = linspace(0,100,10000);
            end
            
            if nargin < 2
                X = zeros(obj.totalPeakParams,1);
            end
            
            %Check that there is a model to fit
            if obj.numPeaks<1
                error('No peaks in model');
            end
            
            %Create a matrix to store the values of each peak
            components = zeros(size(w,1), obj.numPeaks, size(w,2));
            
            %Compute the components for each of the peaks
            for p = 1:obj.numPeaks
                if combo(p)
                    %Get the values of X for a given peak
                    peak_param_values = X(obj.getPeakIdxs(p),:);
                    
                    components(:,p,:) = obj.peakModels(p).modelFun(w,peak_param_values);
                end
            end
            
            %Sum to get full H
            H=squeeze(sum(components,2));
            
            if nargout==2
                components=squeeze(components);
            end
        end
        
        %-----------------------------------
        % MODEL FUNCTION DERIVATIVE (H')
        %-----------------------------------
        function Hx = getHx(obj, X, w, combo)
            %Compute the partial derivative of multipeak function with
            %respect to the state variables in X
            %
            %   Hx = getHx(X, w, combo)
            %
            %   Input:
            %       X: 1 x totalPeakParams double, state vector
            %       w: 1 x W double, frequency vector
            %       combo: 1 x numPeaks logical, selector for active peaks
            %
            %   Output:
            %       Hx: X x w double, multipeak function partial derivatives value
            
            %Set default to all combos being true
            
            %Set default to all combos being true
            if nargin<4 || isempty(combo)
                combo=true(1,obj.totalPeakParams);
            end
            
            %Check that there is a model to fit
            if obj.numPeaks<1
                error('No peaks in model');
            end
            
            %Allocate derivative matrix
            Hx = zeros(length(w), obj.totalPeakParams);
            
            %Loop through peaks and get partial derivatives
            for peak_num = 1:obj.numPeaks
                if combo(peak_num)
                    select_peak = zeros(size(combo));
                    select_peak(peak_num) = 1;
                    
                    H_peak = obj.getH(X, w, select_peak);
                    state_idxs = find(obj.getPeakIdxs(peak_num));
                    
                    %Get the values of X for a given peak
                    peak_param_values = X(state_idxs,:);
                    
                    %Get the current peak model
                    peak_model = obj.peakModels;
                    
                    %For each peak parameter compute the partial derivative
                    for param_num = 1:length(peak_param_values)
                        %Get peak parameter function handle
                        deriv_fun = peak_model(peak_num).modelDerivFuns{param_num};
                        
                        Hx(:,state_idxs(param_num)) = deriv_fun(w(:), peak_param_values, H_peak);
                    end
                end
            end
        end
        
        %-----------------------------------
        % MODEL FUNCTION DERIVATIVE (H')
        %-----------------------------------
        function Hxx = getHxx(obj, X, w, combo)
            %Compute the partial derivative of multipeak function with
            %respect to the state variables in X
            %
            %   Hx = getHx(X, w, combo)
            %
            %   Input:
            %       X: 1 x totalPeakParams double, state vector
            %       w: 1 x W double, frequency vector
            %       combo: 1 x numPeaks logical, selector for active peaks
            %
            %   Output:
            %       Hx: X x w double, multipeak function partial derivatives value
            
            %Set default to all combos being true
            
            %Set default to all combos being true
            if nargin<4 || isempty(combo)
                combo=true(1,obj.totalPeakParams);
            end
            
            %Check that there is a model to fit
            if obj.numPeaks<1
                error('No peaks in model');
            end
            
            %Allocate derivative matrix
            Hxx = zeros(length(w), obj.totalPeakParams, obj.totalPeakParams);
            
            %Loop through peaks and get partial derivatives
            for peak_num = 1:obj.numPeaks
                if combo(peak_num)
                    select_peak = zeros(size(combo));
                    select_peak(peak_num) = 1;
                    
                    H_peak = obj.getH(X, w, select_peak);
                    state_idxs = find(obj.getPeakIdxs(peak_num));
                    
                    %Get the values of X for a given peak
                    peak_param_values = X(state_idxs,:);
                    
                    %Get the current peak model
                    peak_model = obj.peakModels;
                    
                    %For each peak parameter compute the partial derivative
                    for param_num1 = 1:length(peak_param_values)
                        for param_num2 = param_num1:length(peak_param_values)
                            %Get peak parameter function handle
                            deriv_fun = peak_model(peak_num).model2ndDerivFuns{param_num1, param_num2};
                            
                            Hxx(:,state_idxs(param_num1),state_idxs(param_num2)) = deriv_fun(w(:), peak_param_values, H_peak);
                            if param_num1 ~= param_num2
                                Hxx(:,state_idxs(param_num2),state_idxs(param_num1)) = deriv_fun(w(:), peak_param_values, H_peak);
                            end
                        end
                    end
                end
            end
        end
        
        %-----------------------------------
        % MODEL FUNCTION (H) HANDLE
        %-----------------------------------
        function H_handle = getHHandle(obj, ~, ~, ~)
            
            idxs = false(obj.numPeaks, obj.totalPeakParams);
            
            for p = 1:obj.numPeaks
                idxs(p,:)=logical(obj.getPeakIdxs(p));
            end
            
            fhandles = {obj.peakModels.modelFun};
            
            H_handle = @(X,w,combo)sum(combo.*cell2mat(cellfun(@(peakfun, peaknum)peakfun(w',X(idxs(peaknum,:),:)),fhandles, num2cell(1:p), 'UniformOutput',false)),2);
        end
        
        
        
        function H = getH2(obj, X, w, combo)
            
            idxs = false(obj.numPeaks, obj.totalPeakParams);
            
            for p = 1:obj.numPeaks
                idxs(p,:)=logical(obj.getPeakIdxs(p));
            end
            
            fhandles = {obj.peakModels.modelFun};
            
            funvals = cellfun(@(peakfun, peaknum)peakfun(w,X(idxs(peaknum,:),:)),fhandles, num2cell(1:p), 'UniformOutput',false);
            H = sum(reshape(combo,1,1,obj.numPeaks).*cat(3,funvals{:}),3);
        end
        
        %---------------------------------------------
        % MODEL FUNCTION NUMERICAL DERIVATIVE (H')
        %---------------------------------------------
        function Hx = getHxNumerical(obj, X, w, combo, dx)
            %Numerically estimates the partial derivative of multipeak function with
            %respect to the state variables in X
            %
            %   Hx = getHx(X, w, combo, dx)
            %
            %   Input:
            %       X: 1 x totalPeakParams double, state vector
            %       w: 1 x W double, frequency vector
            %       combo: 1 x numPeaks logical, selector for active peaks
            %       dx: double, step size for numerical derivative
            %
            %   Output:
            %       Hx: X x w double, multipeak function partial derivatives value
            
            if nargin<5
                dx = 1e-7;
            end
            
            %Set default to all combos being true
            if nargin<4 || isempty(combo)
                combo=true(1,obj.totalPeakParams);
            end
            
            %Check that there is a model to fit
            if obj.numPeaks<1
                error('No peaks in model');
            end
            
            %Allocate derivative matrix
            Hx = zeros(length(w), obj.totalPeakParams);
            
                       
            %Loop through each parameter and compute Hx
            for pp = 1:obj.totalPeakParams
                Xd = X;
                Xd(pp) = Xd(pp) + dx;
                
                H = obj.getH(X, w, combo);
                Hdx = obj.getH(Xd, w, combo);
                
                Hx(:, pp) = (Hdx - H)/ dx;
            end
        end
        
        %---------------------------------------------
        % PLOT FUNCTION FOR H
        %---------------------------------------------
        function plotObj = plot(obj, X, w, combo, varargin)
            %Plot function for H
            %
            %  plot(obj, X, w, combo, <plot function arguments>)
            %  plot(obj, X, w, peak_num, <plot function arguments>)
            %
            %   Input:
            %       X: 1 x totalPeakParams double, state vector (default: zeros)
            %       w: 1 x W double, frequency vector (default: 1-100 Hz)
            %       combo: 1 x numPeaks logical OR integer, selector for active peaks (default: all peaks)
            
            if nargin < 3 || isempty(w)
                w = linspace(0,100,10000);
            end
            
            if nargin < 4 || isempty(X)
                X = zeros(obj.totalPeakParams, 1);
            else
                if size(X,1) ~= obj.totalPeakParams
                    error('Peak combo size not valid');
                end
            end
            
            if nargin < 2
                combo = true(size(X));
                
            end
            
            %Allow for the peak number to be used as the input
            if isnumeric(combo) && numel(combo)==1
                if combo > obj.numPeaks
                    error('Peak number not valid');
                end
                
                %Create the combo with the selected peaks on
                c = combo;
                combo = false(size(X));
                combo(c, :) = true;
            else
                if isempty(combo)
                    combo = true(size(X));
                else
                    if length(combo) ~= obj.numPeaks
                        error('Peak combo size not valid');
                    end
                end
            end
            
            H = getH(obj, X, w, combo);
            
            plotObj = plot(w, H, varargin{:});
        end
        
        %---------------------------------------------
        % PLOT FUNCTION FOR Hx
        %---------------------------------------------
        function plotObj = plotHx(obj, X, w, combo, varargin)
            %Plot function for Hx
            %
            %  plotHx(obj, X, w, combo, <plot function arguments>)
            %  plotHx(obj, X, w, peak_num, parameter_num, <plot function arguments>)
            %
            %   Input:
            %       X: 1 x totalPeakParams double, state vector (default: zeros)
            %       w: 1 x W double, frequency vector (default: 1-100 Hz)
            %       combo: 1 x numPeaks logical selector for active peaks (default: all peaks)
            
            if nargin < 3 || isempty(w)
                w = linspace(0,100,10000);
            end
            
            if nargin < 4 || isempty(X)
                X = zeros(obj.totalPeakParams, 1);
            else
                if size(X,1) ~= obj.totalPeakParams
                    error('Peak combo size not valid');
                end
            end
            
            %Set default combo to all peaks
            if nargin < 2
                combo = true(size(X));
            end
            
            %Allow for the peak number to be used as the input
            if isnumeric(combo) && numel(combo)==1
                peak_num = combo;
                
                if peak_num > obj.numPeaks
                    error('Peak number not valid');
                end
                
                %Get all the partial derivatives for that peak
                Hx = getHx(obj, X, w, []);
                peakHx = Hx(:, obj.getPeakIdxs(peak_num));
                
                %Check to see if the next argument is a parameter
                if ~isempty(varargin)
                    %Check for it being a string (i.e. plot option)
                    if ~ischar(varargin{1})
                        num_params = obj.peakModels(peak_num).numPeakParams;
                        param_num = varargin{1};
                        
                        %Check if it's a valid parameter number
                        if param_num > num_params
                            error('Invalid peak parameter number');
                        end
                        
                        %Default to all peaks
                        if isempty(param_num)
                            param_num = true(1, num_params);
                        end
                        
                        %Grab the rest of the arguments to pass the plot
                        %function
                        args_rest = {varargin{2:end}};
                        
                        plotObj = plot(w, peakHx(:, param_num), args_rest{:});
                    else
                        plotObj = plot(w, peakHx, varargin{:});
                    end
                    
                end
            else
                %Set default combo to all peaks
                if isempty(combo)
                    combo = true(size(X));
                else
                    if length(combo) ~= obj.numPeaks
                        error('Peak combo size not valid');
                    end
                end
                
                %Plot
                Hx = getHx(obj, X, w, combo);
                plotObj = plot(w, Hx, varargin{:});
            end
        end
    end
end
        