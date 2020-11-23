classdef PeakObjParam < handle
    %PeakObjParam State param object for PeakModel
    %
    %   obj = PeakObjParam(paramName, linkFunClass, linkFunParams)
    %
    %   Input:
    %       paramName: string, parameter name
    %       linkFunClass: string, 'sigmoid' (default), 'exp', 'poly', 'identity'
    %           'smoothdoubledecay', 'linear'
    %       linkFunParams: 1 x P double vector of param values
    %
    %   LINK FUNCTIONS:
    %        Sigmoid ('sigmoid')
    %             Function:
    %                 y(x) = (xmax - xmin) * exp(x)/(1 + exp(x)) + xmin
    %             linkFunParams:
    %                 [xmin, xmax]
    %
    %         Exponential ('exp');
    %             Function:
    %                 y(x) = sign(dir) * exp(x) + offset;
    %             linkFunParams:
    %                 [offset, dir]
    %
    %         Polynomial ('poly');
    %             Function:
    %                 y(x) = sign(dir) * x^p + offset;
    %             linkFunParams:
    %                 [p, offset, dir]
    %
    %          Identity Decay ('identity');
    %             Function:
    %                 y(x) = x;
    %             linkFunParams:
    %                 none
    %
    %   Copyright 2019 Prerau Laboratory
    %
    %   Last modified 5/17/2019 by Michael Prerau
    
    properties (Access = public)
        paramName; %Parameter name
        linkFunClass; %Function class for the link function
        
        linkFun; %Handle to link function
        linkDerivFun; %Handle to link function derivative
        linkDeriv2Fun; %Handle to link function second derivative
        
        linkFunParams; %Link function parameter
        
        linkInvFun; % Handle to inverse of link function
    end
    
    methods (Access = public)
        %-----------------------------------
        %       CONSTRUCTOR
        %-----------------------------------
        function obj = PeakObjParam(paramName, linkFunClass, linkFunParams)
            %Constructor for the PeakObjParam class
            %
            %   obj = PeakObjParam(paramName, linkFunClass, linkFunParams)
            %
            %   Input:
            %       paramName: string, parameter name
            %       linkFunClass: string, 'sigmoid' (default), 'exp', 'identity', 'poly'
            %           'smoothdoubledecay', 'linear'
            %       linkFunParams: 1 x P double vector of param values
            if nargin == 0
                paramName = 'Parameter';
            end
            
            if nargin < 2
                linkFunClass = 'Identity';
            end
            
            if nargin < 3
                linkFunParams = [];
            end
            
            obj.paramName = paramName;
            
            switch lower(linkFunClass)
                %Sigmoid link function
                case {'bound','bounded','sigmoid','sigmoidal'}
                    obj.linkFunClass = 'Sigmoid';
                    
                    %Handle parameters
                    if isempty(linkFunParams)
                        xmin = 0;
                        xmax = 1;
                        
                        linkFunParams = [xmin, xmax];
                    else
                        if length(linkFunParams) ~= 2
                            error('Sigmoid function must have two parameters: max and min');
                        end
                        
                        if linkFunParams(1) >= linkFunParams(2)
                            error('Max value is less  or equal than min value for bound function');
                        end
                    end
                    
                    obj.linkFun = @obj.sigmoidFun;
                    obj.linkDerivFun = @obj.sigmoidFun_deriv;
                    obj.linkDeriv2Fun = @obj.sigmoidFun_deriv2;
                    
                    obj.linkFunParams = linkFunParams;
                    
                    obj.linkInvFun = @obj.sigmoidFun_inv;
                    
                    %Exponential one-sided link function
                case {'one_sided','onesided', 'one sided', 'exp','exponential'}
                    obj.linkFunClass = 'Exponential';
                    
                    %Handle parameters
                    if isempty(linkFunParams)
                        offset=0;
                        dir=1;
                        
                        linkFunParams = [offset, dir];
                    else
                        dir=linkFunParams(2);
                        
                        
                        if length(linkFunParams) ~= 2
                            error('One-sided function must have two parameters: offset and dir');
                        end
                        
                        if abs(dir)~=1
                            error('Direction parameter of exponential one-sided function must be 1 or -1');
                        end
                    end
                    
                    obj.linkFun = @obj.expFun;
                    obj.linkDerivFun = @obj.expFun_deriv;
                    obj.linkDeriv2Fun = @obj.expFun_deriv2;
                    
                    obj.linkFunParams = linkFunParams;
                    
                    obj.linkInvFun = @obj.expFun_inv;
                    
                    %Polynomial link function
                case {'polynomial','poly'}
                    obj.linkFunClass = 'Polynomial';
                    
                    %Handle parameters
                    if isempty(linkFunParams)
                        p = 2;
                        offset = 0;
                        dir = 1;
                        
                        linkFunParams = [p, offset, dir];
                    else
                        dir=linkFunParams(3);
                        
                        
                        if length(linkFunParams) ~= 3
                            error('One-sided function must have three parameters: p, offset, and dir');
                        end
                        
                        if abs(dir)~=1
                            error('Direction parameter of exponential one-sided function must be 1 or -1');
                        end
                    end
                    
                    obj.linkFun = @obj.polynomialFun;
                    obj.linkDerivFun = @obj.polynomialFun_deriv;
                    obj.linkDeriv2Fun = @obj.polynomialFun_deriv2;
                    
                    obj.linkFunParams = linkFunParams;
                    
                    obj.linkInvFun = @obj.polynomialFun_inv;
                    
                    %Identity function
                otherwise
                    obj.linkFunClass = 'Identity';
                    
                    obj.linkFun = @obj.identityFun;
                    obj.linkDerivFun = @obj.identityFun_deriv;
                    obj.linkDeriv2Fun = @obj.identityFun_deriv2;
                    
                    obj.linkFunParams = linkFunParams;
                    
                    obj.linkInvFun = @obj.identityFun_inv;
                    
            end
        end
    end
    
    methods (Access = protected)
        %-----------------------------------
        %  SIGMOID BOUNDING FUNCTION
        %-----------------------------------
        function out = sigmoidFun(obj, x_bar)
            xmin = obj.linkFunParams(1);
            xmax = obj.linkFunParams(2);
            
            out = (xmax-xmin)./(1+exp(-x_bar))+xmin;
        end
        
        %-----------------------------------
        %  SIGMOID BOUNDING DERIVATIVE
        %-----------------------------------
        function out = sigmoidFun_deriv(obj, x_bar)
            xmin = obj.linkFunParams(1);
            xmax = obj.linkFunParams(2);
            
            out = (xmax-xmin)./(1+exp(-x_bar)) .* (1-1./(1+exp(-x_bar)));
        end
        
        %-----------------------------------------
        %  SIGMOID BOUNDING 2nd DERIVATIVE
        %-----------------------------------------
        function out = sigmoidFun_deriv2(obj, x_bar)
            xmin = obj.linkFunParams(1);
            xmax = obj.linkFunParams(2);
            
            sigma = 1./(1+exp(-x_bar));
            out = (xmax-xmin).*(sigma-3*sigma.^2+2*sigma.^3);
        end
        
        %-----------------------------------
        %  SIGMOID BOUNDING INVERSE
        %-----------------------------------
        function out = sigmoidFun_inv(obj, x)
            xmin = obj.linkFunParams(1);
            xmax = obj.linkFunParams(2);
            
            out = -log((xmax-xmin)./(x - xmin)-1);
        end
        
        %-----------------------------------
        %  ONE-SIDED EXPONENTIAL FUNCTION
        %-----------------------------------
        function out = expFun(obj, x_bar)
            offset = obj.linkFunParams(1);
            dir = obj.linkFunParams(2);
            
            out = sign(dir).*exp(x_bar)+offset;
        end
        
        %-----------------------------------
        %  ONE-SIDED EXPONENTIAL DERIVATIVE
        %-----------------------------------
        function out = expFun_deriv(obj, x_bar)
            dir = obj.linkFunParams(2);
            
            out = sign(dir).*exp(x_bar);
        end
        
        %----------------------------------------
        %  ONE-SIDED EXPONENTIAL 2nd DERIVATIVE
        %----------------------------------------
        function out = expFun_deriv2(obj, x_bar)
            dir = obj.linkFunParams(2);
            
            out = sign(dir).*exp(x_bar);
        end
        
        %-----------------------------------
        %  ONE-SIDED EXPONENTIAL Inverse
        %-----------------------------------
        function out = expFun_inv(obj, x)
            offset = obj.linkFunParams(1);
            dir = obj.linkFunParams(2);
            
            out = log(sign(dir).*(x - offset));
        end
        
        %-----------------------------------
        %  POLYNOMIAL FUNCTION
        %-----------------------------------
        function out = polynomialFun(obj, x_bar)
            p = obj.linkFunParams(1);
            offset = obj.linkFunParams(2);
            dir = obj.linkFunParams(3);
            
            out = sign(dir) .* x_bar.^p + offset;
        end
        
        %-----------------------------------
        %  POLYNOMIAL DERIVATIVE
        %-----------------------------------
        function out = polynomialFun_deriv(obj, x_bar)
            p = obj.linkFunParams(1);
            dir = obj.linkFunParams(3);
            
            out = sign(dir).*p.*(x_bar.^(p-1));
        end
        
        %-----------------------------------
        %  POLYNOMIAL 2nd DERIVATIVE
        %-----------------------------------
        function out = polynomialFun_deriv2(obj, x_bar)
            p = obj.linkFunParams(1);
            dir = obj.linkFunParams(3);
            
            out = sign(dir).*p.*(p-1).*(x_bar.^(p-2));
        end
        
        %-----------------------------------
        %  POLYNOMIAL FUNCTION
        %-----------------------------------
        function out = polynomialFun_inv(obj, x)
            p = obj.linkFunParams(1);
            offset = obj.linkFunParams(2);
            dir = obj.linkFunParams(3);
            
            out = nthroot(sign(dir) .* (x - offset),p);
        end
        
        %-----------------------------------
        % IDENTITY FUNCTION
        %-----------------------------------
        function out = identityFun(~, x_bar)
            out = x_bar;
        end
        
        %-----------------------------------
        %  IDENTITY DERIVATIVE
        %-----------------------------------
        function out = identityFun_deriv(~, x_bar)
            out = eye(size(x_bar));
        end
        
        %-----------------------------------
        %  IDENTITY 2nd DERIVATIVE
        %-----------------------------------
        function out = identityFun_deriv2(~, x_bar)
            % out = eye(size(x_bar));
            out = zeros(size(x_bar));
        end
        
        %-----------------------------------
        % IDENTITY INVERSE
        %-----------------------------------
        function out = identityFun_inv(~, x)
            out = x;
        end
        
    end
    
end

