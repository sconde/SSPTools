classdef EmbeddedRK <  SSPTools.Steppers.ERK
    
    properties %(Access = private)
        
        bhat; % embedding weight vector
        isEmbedded = true; % using adaptive step
        isVariableStep = false;
        rejectedStep = 0;
        acceptedStep = 0;
        absTol = 1e-5; % absolute tolerance (used in Embedded-RK);
        relTol = 1.e-4; % relative tolerance (used in Embedded-RK)
        incrFac = 1.2;
        decrFac = 0.2;
        facMax = 0.5;
        facMin = 0.0001;
        phat; % embedding order
        minStepSize;
        maxStepSize;
        safetyFactor;
        nextDt;
        initDt;
        dt_ = NaN;
        t_ = 0;
        badDt_;
        lte_ = 0;
        rejectedDt_ = NaN;
        rejectedLastStep_ = false;
        preditedStepSize_ = NaN; % need to fix this
    end
    
    methods
        function obj = EmbeddedRK(varargin)
            obj = obj@SSPTools.Steppers.ERK(varargin{:});
            
            inpPar = inputParser;
            inpPar.KeepUnmatched = true;
            addParameter(inpPar,'name','Embedded-ERK');
            addParameter(inpPar, 'phat', []); % embedding order
            addParameter(inpPar, 'bhat', []);
            addParameter(inpPar, 'RelTol', 1e-4);
            addParameter(inpPar, 'TolAbs', 1e-5);
            addParameter(inpPar, 'MinStepSize',1e-4);
            addParameter(inpPar, 'MaxStepSize',1e-1);
            addParameter(inpPar, 'StepSizeIncreaseFactor',1.5);
            addParameter(inpPar, 'StepSizeDecreaseFactor',0.5);
            addParameter(inpPar, 'SafetyFactor', 0.8);
            addParameter(inpPar, 'InitialStepSize',1e-4);
            addParameter(inpPar, 'VariableStepSize', false);
            
            inpPar.parse(varargin{:});
            
            obj.bhat = inpPar.Results.bhat;
            obj.phat = inpPar.Results.phat;
            obj.absTol = inpPar.Results.TolAbs;
            obj.relTol = inpPar.Results.RelTol;
            obj.isVariableStep = inpPar.Results.VariableStepSize;
            obj.minStepSize = inpPar.Results.MinStepSize;
            obj.maxStepSize = inpPar.Results.MaxStepSize;
            obj.safetyFactor = inpPar.Results.SafetyFactor;
            obj.incrFac = inpPar.Results.StepSizeIncreaseFactor;
            obj.decrFac = inpPar.Results.StepSizeDecreaseFactor;
            obj.name = inpPar.Results.name;
        end
        
        function [y] = takeStep(obj, dt)
            u0 = obj.u0;
            obj.Y(:,1) = u0;
            
            % intermediate stage value
            for i = 2:obj.s
                temp = u0;
                for j = 1:i-1
                    %keyboard
                    temp = temp + dt*obj.A(i,j)*obj.L(dt + obj.c(j), obj.Y(:,j));
                    %keyboard
                end
                obj.Y(:,i) = temp;
            end
            
            % combine
            y = u0;
            for i = 1:obj.s
                y = y + dt*obj.b(i)*obj.L(dt + obj.c(i), obj.Y(:,i));
            end
            
            % get the embedding
            yhat = u0;
            for i = 1:obj.s
                yhat = yhat + dt*obj.bhat(i)*obj.L(dt + obj.c(i), obj.Y(:,i));
            end
            
            nxDt = stepSizeControl(obj, dt,  y, yhat);
            obj.preditedStepSize_ = nxDt;
        end
        
        
        function dt = startingStepSize(obj)
            % Solving Ordinary Differential Equation I - nonstiff
            % Hairer. pg. 169. Starting Step Size Algorithm
            
            err_fun = @(x, sc_i) sqrt(  sum(((x./sc_i).^2))/length(x)  );
            % a.)
            k1 = obj.L(obj.t, obj.u0);
            sci = obj.absTol + abs(obj.u0)*obj.relTol;
            d0 = err_fun(obj.u0, sci);
            d1 = err_fun(k1, sci);
            
            % b.) first guess for the step size
            if (d0 < 1e-6 || d1 < 1e-6)
                h0 =1e-6;
            else
                h0 = 0.01*(d0/d1);
            end
            
            % c.)
            y1 = obj.u0 + h0*obj.L(obj.t+h0, obj.u0);
            k2 = obj.L(obj.t + h0, y1);
            
            % d.) an estimate of the second derivative
            d2 = (k2 - k1);
            d2 = err_fun(d2, sci)/h0;
            
            % e. )
            max_d1_d2 = max(d1, d2);
            if max_d1_d2 <= 1e-5
                h1 = max(1e-6, 1e-3*h0);
            else
                h1 = (0.01/max_d1_d2)^(1/(obj.p+1));
            end
            
            dt = min(100*h0, h1);
            
        end
        
        function [t, y, dt, err, bad_dt] = getState(obj)
            y = obj.u0;
            t = obj.t_;
            err = obj.lte_;
            dt = obj.preditedStepSize_;
            bad_dt = obj.badDt_;
            
            if ~isempty(obj.dfdx) && (obj.dfdx.systemSize > 1)
                y = reshape(y, obj.dfdx.nx, obj.dfdx.systemSize);
            end
        end
    end
    
    methods ( Access = protected )
        
        function nextDt = stepSizeControl(obj, dt,  y, yhat)
            % Automatic Step Size Control
            % Hairer. Solving ODE I. pg. 167
            
            
            lte = (y - yhat);
            lte = norm(lte, Inf);
            obj.lte_ = lte;
            sci = obj.absTol + max(abs(obj.u0), abs(y))*obj.relTol;
            err = sqrt(sum((lte./sci).^2)/length(y));
            err = lte;
            
            %err = norm(abs(lte), Inf);
            q = min(obj.p, obj.phat);
            
            if obj.rejectedLastStep_
                fac = 1;
            else
                %fac = obj.safetyFactor;
                %fac = (0.25)^(1/(q+1));    % 119 accepted + 14 rejected
                %fac = 0.9;                 % 103 accepted + 75 rejected
                %fac = 0.8;                 % 110 accepted + 26 rejected
                fac = (0.38)^(1/(q+1));     % 110 accepted + 17 rejected
            end
            
            if abs(err) <= 1 % accept the solution
                obj.badDt_ = NaN;
                obj.rejectedLastStep_ = false;
                obj.acceptedStep = obj.acceptedStep + 1;
                obj.t_ = obj.t_ + dt;
                h_new = dt*min(obj.incrFac, fac*(1/err)^(1/(q+1)));
                h_new = dt*obj.incrFac;
                obj.u0 = y; % accept the solution
            else % reject the solution
                obj.badDt_ = dt;
                obj.t_ = obj.t_;
                obj.rejectedLastStep_ = true;
                %h_new = dt*min(obj.incrFac, max(fac*(1/err)^(1/(q+1)),1e-6));
                h_new = dt*min(obj.incrFac, fac*(1/err)^(1/(q)));
                h_new = dt*obj.decrFac;
                obj.rejectedStep = obj.rejectedStep +1;
                obj.u0;
            end
            
            nextDt = h_new;
            obj.lte_ = lte;
        end
    end
end
