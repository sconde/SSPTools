classdef EmbeddedERK < SSPTools.Steppers.ERK
    
    properties
        bhat; % embedding weight vector
        isEmbedded = false; % using adaptive step
        variableStep = false;
        aTol; % absolute tolerance (used in Embedded-RK);
        rTol; % relative tolerance (used in Embedded-RK)
        incrFac = 1.2;
        decrFac = 0.2;
        facMax = 0.5;
        facMin = 0.0001;
        phat; % embedding order
        minStepSize;
        maxStepSize;
        safetyFactor;
        initDt;
        lte_ = 0;
        rejectedDt_ = NaN;
        rejectedLastStep_ = false;
    end
    
    properties ( Access = private)
        isExplicit = true;
        isMSRK = true; 	% Multi-Stage Runge-Kutta
        isButcher = true; % starting with Butcher formulation
        isLowStorage = false;  % need a way to determine is low-storage
        n;
        Y;
    end
    
    methods
        
        function obj = EmbeddedERK(varargin)
            obj = obj@SSPTools.Steppers.ERK(varargin{:});
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'name','Embedded-ERK');
            addParameter(p, 'isSSP', false);
            addParameter(p, 'isButcher', true);
            addParameter(p, 'phat', []); % embedding order
            addParameter(p, 'bhat', []);
            addParameter(p, 'RelTol', 1e-4);
            addParameter(p, 'AbsTol', 1e-5);
            addParameter(p, 'MinStepSize',1e-4);
            addParameter(p, 'MaxStepSize',1e-1);
            addParameter(p, 'StepSizeIncreaseFactor',1.2);
            addParameter(p, 'StepSizeDecreaseFactor',0.5);
            addParameter(p, 'SafetyFactor', 0.8);
            addParameter(p, 'InitialStepSize',1e-4);
            addParameter(p, 'VariableStepSize', false);
            p.parse(varargin{:});
            
            obj.variableStep = p.Results.VariableStepSize;
            keyboard
                        
            % embedded-rk parameters
            if ~isempty(p.Results.bhat)
                obj.bhat = p.Results.bhat;
                obj.isEmbedded = true;
                obj.phat = p.Results.phat;
                obj.minStepSize = p.Results.MinStepSize;
                obj.maxStepSize = p.Results.MaxStepSize;
                obj.incrFac = p.Results.StepSizeIncreaseFactor;
                obj.decrFac = p.Results.StepSizeDecreaseFactor;
                obj.safetyFactor = p.Results.SafetyFactor;
                obj.initDt = p.Results.InitialStepSize;
                obj.dt_ = obj.initDt;
                %obj.nextDt = -Inf;
            end
            
            obj.isEmbedded = obj.isEmbedded && obj.variableStep;
            
            obj.rTol = p.Results.RelTol;
            obj.aTol = p.Results.AbsTol; % why is this give me the butcher matrix?
            
            keyboard
            
        end % end constructor
        
        
    end
    
    methods %( Access = protected )
        
        function [y, dt] = takeStep(obj, dt)
            % function [y, dt] = takeStep(dt)
            % returns the new solution (y) and time-step taken (dt)
            u0 = obj.u0;
            obj.Y(:,1) = u0;
            obj.dt_ = dt;
            
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
            
            % compute the embedding solution
            if obj.isEmbedded % just to make sure
                yhat = u0;
                for i = 1:obj.s
                    yhat = yhat + dt*obj.bhat(i)*obj.L(dt + obj.c(i), obj.Y(:,i));
                end
                
                [y, t] = obj.stepSizeControl(dt, y, yhat);
            else
                obj.nextDt = dt;
                t = obj.t + dt;
            end
            
            obj.u0 = y;
            obj.t  = t;
        end
        
        function [t, y, nextDt, err, rejectedDt] = getState(obj)
            y = obj.u0;
            t = obj.t;
            err = obj.lte_;
            rejectedDt = obj.rejectedDt_;
            %keyboard
            if (~obj.isEmbedded) || (t == 0)
                nextDt = obj.dt_;
%                 if ~isempty(obj.dt_) % use the user starting step size
%                     nextDt = obj.dt_;
%                 else % use the estimated starting step size
%                     nextDt = obj.startingStepSize();
%                 end
            else
                nextDt = obj.nextDt;
            end
                        
            if ~isempty(obj.dfdx) && (obj.dfdx.systemSize > 1)
                y = reshape(y, obj.dfdx.nx, obj.dfdx.systemSize);
            end
        end
        
        function dt = startingStepSize(obj)
            % Solving Ordinary Differential Equation I - nonstiff
            % Hairer. pg. 169. Starting Step Size Algorithm
            
            err_fun = @(x, sc_i) sqrt(  sum(((x./sc_i).^2))/length(x)  );
            % a.)
            k1 = obj.L(obj.t, obj.u0);
            sci = obj.aTol + abs(obj.u0)*obj.rTol;
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
        
    end
    
    methods ( Access = protected )
        
        
        function [y, t] = stepSizeControl(obj,dt,  y, yhat)
            % Automatic Step Size Control
            % Hairer. Solving ODE I. pg. 167
            
            keyboard
            lte = (y - yhat);
            sc_i = obj.aTol + max(abs(obj.u0), abs(y))*obj.rTol;
            %sc_i = obj.aTol + max(abs(yhat), abs(y))*obj.rTol;
            %err = sqrt(sum((lte./sc_i).^2)/length(lte));
            err = norm(abs(lte), Inf);
            
            if err < 1
                % accept the solution
                q = min(obj.p, obj.phat);
                dt_op = dt*(1/err).^(1/q);
                t = obj.t + dt;
            else
                y = obj.u0;
                t = obj.t;
            end
            
            dt_new = dt*min(obj.incrFac, max(obj.decrFac, ...
                obj.safetyFactor*dt_op));
            dt = dt_new;
            obj.nextDt = dt;
        end
    end
    
    
end
