classdef EmbeddedRK <  SSPTools.Steppers.ERK
    
    properties (Access = protected)%(Access = private)   
        bhat; % embedding weight vector
        isEmbedded = true; % using adaptive step
        isVariableStep = false;
        rejectedStep = 0;
        acceptedStep = 0;
        absTol;% absolute tolerance (used in Embedded-RK);
        relTol;% relative tolerance (used in Embedded-RK)
        incrFac;
        decrFac;
        facMax;
        facMin;
        safe;
        fac1 = 0.2; % facMin
        fac2 = 10.0; % facMax;
        beta_;
        expo1;
        facc1;
        facc2;
        facold = 1.0E-4;
        phat; % embedding order
        minStepSize;
        maxStepSize;
        safetyFactor;
    end
    
    properties(Access = protected)
        initDt;
        t_ = 0;
        badDt_;
        lte_ = 0;
        rejectedDt_ = NaN;
        rejectedLastStep_ = false;
        preditedStepSize_ = NaN; % need to fix this
        fac11;
        last = 0;
        naccpt = 0;
        nrejct = 0;
        reject = 0;
        oldT = 0;
        tFinal;
        dtMax;
        posneg;
        nfcn = 0;
        nstep = 0;
    end
    
    methods
        function obj = EmbeddedRK(varargin)
            obj = obj@SSPTools.Steppers.ERK(varargin{:});
            
            inpPar = inputParser;
            inpPar.KeepUnmatched = true;
            addParameter(inpPar,'name','Embedded-ERK');
            addParameter(inpPar, 'phat', []); % embedding order
            addParameter(inpPar, 'bhat', []);
            addParameter(inpPar, 'RelativeTolerance', 1e-7);
            addParameter(inpPar, 'AbsoluteTolerance', 1e-7);
            addParameter(inpPar, 'MinStepSize',1e-4);
            addParameter(inpPar, 'MaxStepSize',[]);
            addParameter(inpPar, 'StepSizeIncreaseFactor',1.5);
            addParameter(inpPar, 'StepSizeDecreaseFactor',0.5);
            addParameter(inpPar, 'SafetyFactor', 0.8);
            addParameter(inpPar, 'InitialStepSize',1e-4);
            addParameter(inpPar, 'VariableStepSize', false);
            addParameter(inpPar, 'Tfinal', 1);
            addParameter(inpPar, 'Beta', 0.04);
            addParameter(inpPar, 'Safety', 0.9);
           
            inpPar.parse(varargin{:});
            
            obj.bhat = inpPar.Results.bhat;
            obj.phat = inpPar.Results.phat;
            obj.absTol = inpPar.Results.AbsoluteTolerance;
            obj.relTol = inpPar.Results.RelativeTolerance;
            obj.isVariableStep = inpPar.Results.VariableStepSize;
            obj.minStepSize = inpPar.Results.MinStepSize;
            obj.maxStepSize = inpPar.Results.MaxStepSize;
            obj.safetyFactor = inpPar.Results.SafetyFactor;
            obj.incrFac = inpPar.Results.StepSizeIncreaseFactor;
            obj.decrFac = inpPar.Results.StepSizeDecreaseFactor;
            obj.name = inpPar.Results.name;
            obj.beta_ = inpPar.Results.Beta;
            obj.safe = inpPar.Results.Safety;
            
            obj.expo1 = 0.2 - obj.beta_ * 0.75;
            obj.facc1 = 1.0 / obj.fac1 ;
            obj.facc2 = 1.0 / obj.fac2 ;
            obj.tFinal = inpPar.Results.Tfinal;
            
            % set maximum step size allow
            if isempty(obj.maxStepSize)
                obj.dtMax = obj.tFinal - obj.t;
            else
                obj.dtMax = obj.maxStepSize;
            end
            
            % define the posneg variable
            if obj.dtMax < 0
                obj.posneg = -1;
            else
                obj.posneg = 1;
            end
            
        end
        
        function [y] = takeStep(obj, dt)
            % not implemented
            error('Not implemented');
        end
        
        
        function summary(obj)
            fprintf(1, 'rtol= %12.5e, atol= %12.5e, fcn = %d, step = %d, accpt = %d, reject = %d\n',...
                obj.relTol, obj.absTol, obj.nfcn, obj.nstep, obj.naccpt, obj.nrejct);
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
            obj.nextDt = dt;
            obj.nfcn = 2;
        end
        
        function [t, y, dt, err, bad_dt] = getState(obj)
            
            y = obj.u0;
            t = obj.t;
            err = norm(obj.lte_, Inf);
            
            if isempty(obj.nextDt) && isempty(obj.dt_)
                obj.startingStepSize();
            end
            
            dt = obj.nextDt;
            
            if (obj.naccpt && obj.reject)
                bad_dt = obj.dt_;
            else
                bad_dt = NaN;
            end
            
            if ~isempty(obj.dfdx) && (obj.dfdx.systemSize > 1)
                y = reshape(y, obj.dfdx.nx, obj.dfdx.systemSize);
            end
        end
    end
    
    methods ( Access = protected )
        
        function sk = sci_fun(obj, y, yhat)
            sk = obj.absTol + obj.relTol * max(abs(obj.u0), abs(yhat));
        end
        
        function [y, t] = stepSizeControl(obj, dt,  y, yhat)
            % Automatic Step Size Control
            % Hairer. Solving ODE I. pg. 167
            
            
            h = dt;
            lte = abs(y - yhat);
            obj.lte_ = norm(lte, Inf);
            
            sk = obj.sci_fun(y, yhat);
            err = sum((lte./sk).^2);
            err = sqrt(err /length(y));
            
            %/* computation of hnew */
            obj.fac11 = err.^(obj.expo1);
            
            % /* Lund-stabilization */
            fac = obj.fac11/(obj.facold.^(obj.beta_));
            
            % /* we require fac1 <= hnew/h <= fac2 */
            fac = max(obj.facc2, min(obj.facc1, fac/obj.safe));
            hnew = h / fac;
            
            if abs(err) <= 1 % accept the solution
                % accept the step
                obj.facold = max(err, 1.0E-4);
                obj.naccpt = obj.naccpt + 1;
                obj.oldT = obj.t;
                t = obj.t + dt;
                
                % /* normal exit */
                if (obj.last)
                    obj.nextDt = hnew;
                    %obj.t = t;
                    idid = 1;
                    keyboard
                    %break;
                end
                
                if (abs(hnew) > obj.dtMax)
                    hnew = obj.posneg * obj.dtMax;
                end
                
                if (obj.reject)
                    hnew = obj.posneg * min(abs(hnew), abs(h));
                end
                
                obj.reject = 0;
                
            else % reject the solution
                % reject
                
                hnew = h/(min(obj.facc1, obj.fac11/obj.safe));
                obj.reject = 1;
                y = obj.u0;
                t = obj.t;
                
                if (obj.naccpt >= 1)
                    obj.nrejct = obj.nrejct + 1;
                end
                obj.last = 0;
                obj.reject = 1;
                obj.nextDt = hnew;
            end
            
            obj.lte_ = lte;
            obj.nextDt = hnew;
            
        end
    end
end
