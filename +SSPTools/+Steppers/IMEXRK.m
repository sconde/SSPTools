
classdef IMEXRK < SSPTools.Steppers.RK
    
    properties
        At; bt; ct; dgdx; dgdt; ImplicitProblem; rt;
        pim;
    end
    
    properties ( Access = private)
        isExplicit = false;
        isMSRK = true; 	% Multi-Stage Runge-Kutta
        isButcher = true; % starting with Butcher formulation
        isLowStorage = false;  % need a way to determine is low-storage
        n;
        Y;
        F;
        G;
        Fvec;
        Gvec;
        isImplicitLinear = true; %so far handling linear advection
        NL;
        DT;
        solver;
        I;
        utex;
        utimp;
    end
    
    methods
        
        function obj = IMEXRK(varargin)
            obj = obj@SSPTools.Steppers.RK(varargin{:});
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'name','IMEXRK');
            addParameter(p, 'isSSP', false);
            addParameter(p, 'isButcher', true);
            addParameter(p, 'ImplicitODE', []);
            addParameter(p, 'At', []);
            addParameter(p, 'bt', []);
            addParameter(p, 'rt', []);
            addParameter(p, 'pim', []);
            addParameter(p, 'isLowStorage', false);
            addParameter(p, 'dgdx', []);
            addParameter(p, 't', 0.0);
            p.parse(varargin{:});
            
            obj.At = p.Results.At;
            obj.bt = p.Results.bt;
            obj.ct = sum(obj.A,2);
            obj.name = p.Results.name;
                        
            obj.F = zeros(obj.n, obj.s);
            obj.G = zeros(obj.n, obj.s);
            obj.t = p.Results.t;
            
            % get the SSP coefficient of the method
            if ~isempty(p.Results.rt) && ~isinf(p.Results.rt)
                obj.rt = p.Results.rt;
                obj.isSSP = (obj.isSSP && true);
            else
                %TODO: this only get the RK method's SSP, not IMEXRK ssp coefficient.
                % need to calculate this correctly
                obj.rt = obj.am_radius(obj.At, obj.bt(:));
                if obj.rt > 0
                    obj.isSSP = (obj.isSSP && true);
                else
                    obj.isSSP = false;
                end
            end
            
            % check that r and rt are both positive
            
            if ~isempty(p.Results.pim)
                obj.pim = p.Results.pim;
            else
                obj.pim = obj.p;
            end
            
            
            if ~isempty(p.Results.dgdx)
                obj.dgdx = p.Results.dgdx;
                obj.ImplicitProblem = obj.dgdx.problem;
            elseif ~isempty(p.Results.ImplicitODE)
                obj.ImplicitProblem = p.Results.ImplicitODE;
            end
            
            if isa(obj.ImplicitProblem, 'TestProblems.ODEs.ODE') && ...
                    isa(obj.ExplicitProblem, 'TestProblems.ODEs.ODE')
                
                % assume the implicit problem is going to be a nonlinear
                % problem
                obj.solver = @(y, dt, i) nonlinearImplicitStage( obj, y, dt, i );
                obj.F = @(t,y) obj.ExplicitProblem.f(t,y);
                obj.G = @(t,y) obj.ImplicitProblem.f(t,y);
                
            elseif isa(obj.ImplicitProblem, 'TestProblems.PDEs.LinearDiffusion')
                %FIX : this is only working for the diffusion paraemter
                obj.solver = @(y, dt, i) linearSolve(obj, y, dt,i);
                obj.F = @(t,y) obj.dfdx.L(t,y);
                obj.G = @(t,y) obj.dgdx.L(t,y);
                obj.DT = -obj.ImplicitProblem.nu*obj.dgdx.D;
            elseif isa(obj.ImplicitProblem, 'TestProblems.PDEs.LinearAdvection')
                obj.solver = @(y, dt, i) linearSolve(obj, y, dt,i);
                obj.F = @(t,y) obj.dfdx.L(t,y);
                obj.G = @(t,y) obj.dgdx.L(t,y);
                obj.DT = obj.ImplicitProblem.a*obj.dgdx.D;
            else
                obj.solver = @(y, dt, i) nonlinearImplicitStage( obj, y, dt, i );
                obj.F = @(t,y) obj.dfdx.L(t,y);
                obj.G = @(t,y) obj.dgdx.L(t,y);
            end
            
            if isa(obj.y0, 'function_handle')
                obj.u0 = obj.y0(obj.x);
            else
                obj.u0 = obj.y0;
            end
           
            if obj.isSSP
                obj.name = sprintf('SSP%d(%d,%d,%d)%d',...
                   obj.p,obj.s, obj.s,obj.pim, obj.plin);
            else
                obj.name = sprintf('IMEX%d(%d,%d,%d)%d',...
                   obj.p,obj.s, obj.s,obj.pim, obj.plin);
            end
            
            obj.n = size(obj.u0,1);
            obj.I = speye(obj.n);
            obj.Y = zeros(obj.n, obj.s);
            obj.Gvec = zeros(obj.n, obj.s);
            obj.Fvec = zeros(obj.n, obj.s);
        end
        
        
    end
    
    methods
        
        function butcherCoef(obj)
            % print the explicit method
            fprintf('Explicit RK\n');
            obj.printCoeff(obj.A, obj.b, obj.c);
            
            % print the implicit method
            fprintf('\n\nImplicit RK\n');
            obj.printCoeff(obj.At, obj.bt, obj.ct);
        end
        
        function [y, dt] = takeStep(obj, dt)
            % function [y, dt] = takeStep(dt)
            % returns the new solution (y) and time-step taken (dt)
            
            u0 = obj.u0;
            
            % first stage implicit solve
            te = obj.solver(u0,dt, 1);
            obj.Fvec(:, 1) = obj.F(obj.t + dt*obj.c(1), te);
            obj.Gvec(:, 1) = obj.G(obj.t + dt*obj.ct(1), te);
            
            % intermediate stage value
            for i = 2:obj.s
                
                tempt = u0;
                for j = 1:i-1
                    tempt = tempt + dt*obj.A(i,j)*obj.Fvec(:,j) + ...
                        dt*obj.At(i,j)*obj.Gvec(:,j);
                end
                te = obj.solver(tempt, dt, i);
                obj.Fvec(:, i) = obj.F(obj.t + dt * obj.c(i), te);
                obj.Gvec(:, i) = obj.G(obj.t + dt * obj.ct(i), te);
            end

            % combine
            y = u0 + dt*obj.Fvec*obj.b(:) + dt*obj.Gvec*obj.bt(:);
            
            obj.u0 = y;
            obj.t = obj.t + dt;
        end


        function stabilityFunction(obj)
    		error('not yet implemented');    
        end
        
    end
    
    methods (Access = private)
        
        function  y = linearSolve(obj, y, dt,i)
            y = (obj.I - dt*obj.At(i,i)*obj.DT)\y;
            
        end
        
        function y = nonlinearImplicitStage( obj, y, dt, i )
            % not right
            epss = 1e-14;
            options = optimset('Display','off', 'TolFun',epss, 'TolX',epss);
            fzero = @(K) (y + dt*obj.At(i,i)*obj.G(dt + obj.ct(i), K)  - K);
            [Y, fval] = fsolve(@(Y) fzero(Y), y, options);
            y = Y;
        end
        
    end
    
    
end
