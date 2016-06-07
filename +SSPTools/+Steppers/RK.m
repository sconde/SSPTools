classdef RK < handle
    
    properties
        name;			% Name of time-stepping method.
        dfdt;
        dfdx;
        ExplicitProblem;
        A = []; b = []; c = []; alpha = []; s = [];
        r;
        p; % Nonlinear-Order
        plin; % Linear-Order
        t0;
        y0;
        t;
        isSSP;
        bhat; % embedding weight vector
        isEmbedded = false; % using adaptive step
        variableStep = false;
    end
    
    properties (Access = protected)
        steps = 1;
        L;
        CFL;
        dx;
        isLinear;
        haveExact = false;
        x;
        u0;
        systemSize;
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
        nextDt;
        initDt;
        dt_;
    end
    
    methods
        function obj = RK(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'name','MSRK');
            addParameter(p, 'dfdx', []);
            addParameter(p, 'dfdt', []);
            addParameter(p, 'A', []);
            addParameter(p, 'b', []);
            addParameter(p, 'bhat', []);
            addParameter(p, 'r', -Inf);
            addParameter(p, 'p', []);
            addParameter(p, 'phat', []); % embedding order
            addParameter(p, 'plin', []);
            addParameter(p, 'alpha', []);
            addParameter(p, 't', 0.0);
            addParameter(p, 't0', 0);
            addParameter(p, 'y0', []);
            addParameter(p, 'ODE', []);
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
            
            if isa(p.Results.dfdt, 'function_handle')
                obj.dydt = p.Results.dfdt;
            end
            
            if isa(p.Results.y0, 'function_handle')
                obj.y0 = p.Results.y0;
                obj.haveExact = true;
            else
                obj.y0 = p.Results.y0(:);
            end
            
            obj.t = p.Results.t;
            obj.A = p.Results.A;
            obj.b = p.Results.b;
            obj.c = sum(obj.A,2);
            obj.s = numel(obj.b); %infer the number of stages from the length of b
            obj.p = p.Results.p;
            obj.variableStep = p.Results.VariableStepSize;
            
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
            
            % get the SSP coefficient of the method
            if ~isempty(p.Results.r) && ~isinf(p.Results.r)
                obj.r = p.Results.r;
                obj.isSSP = true;
            else
                obj.r = obj.am_radius(obj.A, obj.b(:));
                if obj.r > 0
                    obj.isSSP = true;
                else
                    obj.isSSP = false;
                end
            end
            
            if isempty(p.Results.plin)
                obj.plin = obj.p;
            else
                obj.plin = p.Results.plin;
            end
            
            obj.alpha = p.Results.alpha;
            obj.name = p.Results.name;
            obj.dfdx = p.Results.dfdx;
            obj.t0 = p.Results.t0;
            
            assert(isequal(obj.s, size(p.Results.A,1)),...
                sprintf('RK A:Stage-count -- Num-Rows(A) != %d',obj.s));
            
            if ~isempty(obj.dfdx)
                obj.ExplicitProblem = obj.dfdx.problem;
            elseif ~isempty(p.Results.ODE)
                obj.ExplicitProblem = p.Results.ODE;
            end
            
            if isa(obj.ExplicitProblem, 'TestProblems.ODEs.ODE')
                obj.L = @(t,y) obj.ExplicitProblem.f(t,y);
                obj.CFL = [];
                obj.dx = [];
                obj.x = [];
            else
                obj.L = @(t,y) obj.dfdx.L(t, y);
                obj.dx = obj.dfdx.dx;
                obj.x = p.Results.dfdx.x(:); % make it a vector
            end
            
            obj.isLinear = obj.ExplicitProblem.isLinear;
            
            if obj.isSSP
                obj.name = sprintf('SSP(%d,%d)%d',obj.s, obj.p, obj.plin);
            else
                obj.name = sprintf('RK(%d,%d)%d',obj.s, obj.p, obj.plin);
            end
            
            obj.rTol = p.Results.RelTol;
            obj.aTol = p.Results.AbsTol;
            
        end % RK constructor
    end
    
    methods %( Access = protected )
        
        function butcherCoef(obj)
            %TODO: don't print the zeros
            if obj.isEmbedded
                obj.printCoeff(obj.A, obj.b, obj.c, obj.bhat);
            else
                obj.printCoeff(obj.A, obj.b, obj.c);
            end
        end
        
        function [y, dt] = takeStep(obj, dt) end
        
        function resetInitCondition(obj)
            if isa(obj.y0, 'function_handle')
                obj.u0 = obj.y0(obj.x);
            else
                obj.u0 = obj.y0;
            end
            obj.t = 0.0;
            %obj.nextDt = -Inf;
        end
        
        function [t, y, nextDt] = getState(obj)
            y = obj.u0;
            t = obj.t;
            if (~obj.isEmbedded) || (t == 0)
                nextDt = obj.dt_;
            else
                nextDt = obj.nextDt;
            end
                        
            if ~isempty(obj.dfdx) && (obj.dfdx.systemSize > 1)
                y = reshape(y, obj.dfdx.nx, obj.dfdx.systemSize);
            end
            %             %Q = [r ru E];
            %             density   = u(1:obj.N);
            %             momentum  = u(obj.N+1:2*(obj.N));
            %             energy    = u(2*(obj.N)+1:3*(obj.N));
            %             keyboard
        end
        
        function r = am_radius(obj, A,b)
            %function r = am_radius(A,b)
            %
            %By David Ketcheson
            %
            %Evaluates the Radius of absolute monotonicity
            %of a Runge-Kutta method, given the Butcher array.
            %
            %For an m-stage method, A should be an m x m matrix
            %and b should be a column vector of length m.
            %
            %Accuracy can be changed by modifying the value of eps.
            %Methods with very large radii of a.m. (>50) will require
            %rmax to be increased.
            
            rmax=50; eps=1.e-12;
            
            m=length(b); e=ones(m,1);
            K=[A;b'];
            rlo=-50; rhi=rmax;
            
            while rhi-rlo>eps  %use bisection
                r=0.5*(rhi+rlo);
                X=eye(m)+r*A; beta=K/X; ech=r*K*(X\e);
                if (min(beta(:))<-3.e-16 || max(ech(:))>1.+3.e-16)
                    rhi=r;
                else
                    rlo=r;
                end
            end
            
            if rhi==rmax % r>=rmax
                %error('Error: increase value of rmax in am_radius.m');
                r = Inf;
                %break
            else
                r=rlo;
            end
            
        end
    end
    
    
    methods ( Access = protected )
        
        function [y, t] = stepSizeControl(obj,dt,  y, yhat)
            % Automatic Step Size Control
            % Hairer. Solving ODE I. pg. 167
            
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
        
        function printCoeff(obj, A, b, c, varargin)
            %TODO: don't print the zeros
            del = repmat(' %5.4f ',1, obj.s);
            fprintf(1,['%5.4f |' del '\n'],[c A]');
            fprintf(1,'%s\n',repmat('-',1,8*(obj.s+1)));
            fprintf(1,'%6s |',repmat(' ',1,5));
            fprintf(1,[del '\n'], b);
            if ~isempty(varargin) % method is embedded
                fprintf(1,'%6s |',repmat(' ',1,5));
                fprintf(1,[del '\n'], varargin{1});
            end
        end
        
        function obj = setL(obj)
            obj.L = @(t, y) obj.L(y);
        end
        
    end
end
