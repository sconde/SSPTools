classdef EmbeddedERK < SSPTools.Steppers.EmbeddedRK
    
    properties
        %bhat; % embedding weight vector
    end
    
    properties ( Access = private)
        %isEmbedded = false; % using adaptive step
        variableStep = false;
        aTol = 1e-7; % absolute tolerance (used in Embedded-RK);
        rTol = 1e-7; % relative tolerance (used in Embedded-RK)
        %incrFac = 1.2;
        %decrFac = 0.2;
        %facMax = 0.5;
        %facMin = 0.0001;
        %phat; % embedding order
        %minStepSize;
        %maxStepSize;
        %safetyFactor;
        %initDt;
        %lte_ = 0;
        %rejectedDt_ = NaN;
        %rejectedLastStep_ = false;
        
      
    end
    
    methods
        
        function obj = EmbeddedERK(varargin)
            obj = obj@SSPTools.Steppers.EmbeddedRK(varargin{:});
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'name','Embedded-ERK');
            addParameter(p, 'isSSP', false);
            addParameter(p, 'isButcher', true);
            addParameter(p, 'phat', []); % embedding order
            addParameter(p, 'bhat', []);
            addParameter(p, 'RelTol', 1e-7);
            addParameter(p, 'AbsTol', 1e-7);
            addParameter(p, 'MinStepSize',1e-4);
            addParameter(p, 'MaxStepSize',1e-1);
            addParameter(p, 'StepSizeIncreaseFactor',1.2);
            addParameter(p, 'StepSizeDecreaseFactor',0.5);
            addParameter(p, 'SafetyFactor', 0.8);
            addParameter(p, 'InitialStepSize',1e-4);
            addParameter(p, 'VariableStepSize', false);
            
            p.parse(varargin{:});
            
            obj.variableStep = p.Results.VariableStepSize;
            
            obj.name = p.Results.name;
            assert(~isempty(p.Results.bhat),'Need Bhat Vector');
            
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
            obj.isEmbedded = obj.isEmbedded && obj.variableStep;
            obj.aTol = 1e-7; % absolute tolerance (used in Embedded-RK);
            obj.rTol = 1e-7;

        end % end constructor
        
        
    end
    
    methods %( Access = protected )
        
        
        function [y, dt] = takeStep(obj, dt)
            % function [y, dt] = takeStep(dt)
            % returns the new solution (y) and time-step taken (dt)
            
            u0 = obj.u0;
            obj.Y(:,1) = u0;
            obj.dt_ = dt;
            obj.nstep = obj.nstep + 1;
            
            % intermediate stage value
            for i = 2:obj.s
                temp = u0;
                for j = 1:i-1
                    temp = temp + dt*obj.A(i,j)*obj.L(dt + obj.c(j), obj.Y(:,j));
                end
                obj.Y(:,i) = temp;
            end
            
            obj.nfcn = obj.nfcn + (obj.s-1);
            
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
                
%         function h = startingStepSize(obj) %vdp, y, hmax, atol, rtol,p)
%             % this is good
%             x = obj.t; %starting
%             y = obj.u0;
%             hmax = obj.dtMax;
%             f0 = obj.L(obj.t, obj.u0);
%             % hinit (starting step size)
%             dnf = 0;
%             dny = 0;
%             sk = obj.sci_fun(y, zeros(size(y)));
%             
%             %sk = atol + rtol*abs(y);
%             
%             dnf = sum((f0./sk).^2); % good
%             dny = sum((y./sk).^2); % good
%             
%             % good
%             if ( dnf <= 1e-10 || dny <= 1e-10)
%                 h = 1e-6;
%             else
%                 h = 0.01*sqrt(dny/dnf);
%             end
%             
%             h = min(h, hmax);
%             
%             % perform the explicit euler step
%             yy1 = y + h*f0;
%             f1 = obj.L(x+h, yy1); % fcn (n, x+h, yy1, f1); f1 = vdp(yy1);
%             
%             % /* estimate the second derivative of the solution */
%             %sk = atol + rtol*abs(y);
%             sk = obj.sci_fun(y, zeros(size(y)));
%             sqr = (f1 - f0)./sk;
%             der2  = sum(sqr.^2);
%             der2 = sqrt(der2)/h;
%             
%             % /* step size is computed such that h**iord * max_d(norm(f0),norm(der2)) = 0.01 */
%             der12 = max(abs(der2), sqrt(dnf));
%             
%             if (der12 <= 1e-15)
%                 h1 = max(1e-6, abs(h)*1e-3);
%             else
%                 h1 = (0.01/der12)^(1/max(obj.p, obj.phat));
%             end
%             
%             h = min(100*h, min(h1, hmax));
%             
%         end
        
        function dt = hStepInit(obj)
            % Solving Ordinary Differential Equation I - nonstiff
            % Hairer. pg. 169. Starting Step Size Algorithm
            
            err_fun = @(x, sc_i) sqrt(  sum(((x./sc_i).^2))/length(x)  );
            % a.)
            keyboard
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
        
        
%         function [y, t] = stepSizeControl(obj,dt,  y, yhat)
%             % Automatic Step Size Control
%             % Hairer. Solving ODE I. pg. 167
%             
%             lte = (y - yhat);
%             sk = obj.sci_fun(obj.u0, yhat);
%             err = sum((lte./sk).^2);
%             err = sqrt(err /length(y));
%             
%             keyboard
%             %/* computation of hnew */
%             fac11 = err.^(expo1);
%             
%             % /* Lund-stabilization */
%             fac = fac11/(facold.^(beta));
%             
%             % /* we require fac1 <= hnew/h <= fac2 */
%             fac = max(facc2, min(facc1, fac/safe));
%             hnew = h / fac;
%             
%             if err < 1
%                 % accept the solution
%                 q = min(obj.p, obj.phat);
%                 dt_op = dt*(1/err).^(1/q);
%                 t = obj.t + dt;
%             else
%                 y = obj.u0;
%                 t = obj.t;
%             end
%             
%             dt_new = dt*min(obj.incrFac, max(obj.decrFac, ...
%                 obj.safetyFactor*dt_op));
%             dt = dt_new;
%             obj.nextDt = dt;
%         end
    end
    
    
end
