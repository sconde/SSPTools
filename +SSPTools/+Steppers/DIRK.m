classdef DIRK < SSPTools.Steppers.RK
    
    properties
        
    end
    
    properties ( Access = private)
        isExplicit = false;
        isMSRK = true; 	% Multi-Stage Runge-Kutta
        isButcher = true; % starting with Butcher formulation
        isLowStorage = false;  % need a way to determine is low-storage
        n;
        Y;
        solver;
        NL;
        I;
        isImplicitLinear;
        DT;
    end
    
    methods
        
        function obj = DIRK(varargin)
            obj = obj@SSPTools.Steppers.RK(varargin{:});
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'name','MSRK-DIRK');
            addParameter(p, 'isSSP', false);
            addParameter(p, 'isButcher', true);
            addParameter(p, 'isLowStorage', false);
            addParameter(p, 't', 0.0);
            p.parse(varargin{:});

            obj.name = p.Results.name;  
            obj.t = p.Results.t;
            
            obj.isImplicitLinear = obj.ExplicitProblem.isLinear;
            
            if isa(obj.y0, 'function_handle')
                obj.u0 = obj.y0(obj.x);
                obj.n = size(obj.x,1);
            else
                obj.u0 = obj.y0;
                obj.n = size(obj.u0,1);
            end
            
            obj.name = p.Results.name;            
            obj.Y = zeros(obj.n, obj.s);
            obj.I = speye(obj.n);
            
            if isa(obj.ExplicitProblem, 'TestProblems.ODEs.ODE')
                obj.solver = @(y,dt,i) nonlinearImplicitStage( obj, y, dt, i );
            else
                if obj.isImplicitLinear
                    obj.solver = @(y, dt, i) linearSolve(obj, y, dt,i);
                    obj.DT = obj.dfdx.D;
                    obj.NL = @(t,y) obj.dfdx.L(obj.ExplicitProblem.f(t,y));
                else
                    obj.NL = @nonlinearImplicitStage;
                    obj.solver = @(y, dt, i) nonlinearImplicitStage( obj, y, dt, i );
                end
            end
            
            assert(~isempty(obj.solver),'Implicit Solver is empty');
            
            if obj.isSSP
                obj.name = sprintf('SSP(%d,%d)%d',obj.s, obj.p, obj.plin);
            else
                obj.name = sprintf('RK(%d,%d)%d',obj.s, obj.p, obj.plin);
            end
        end
        
        
    end
    
    methods %( Access = protected )


        function [y] = takeStep(obj, dt)
            
%             %check to see if CFL violation
%             assert((dt/obj.dx) <= obj.CFL, ...
%                 sprintf('ERK: CFL Violation (CFL = %3.2f )',dt/obj.dx) );
            
            u0 = obj.u0;
            
            % first stage implicit solve
            obj.Y(:,1) = obj.solver(u0,dt, 1);
            
            % intermediate stage value
            for i = 1:obj.s
                
                %obj.G(:,i) = obj.dgdx.L(obj.ImplicitProblem.f(dt + obj.ct(i), obj.Y(:,i)));

                temp = u0;
                for j = 1:i-1
                    temp = temp + dt*obj.A(i,j)*obj.L(dt + obj.c(j), obj.Y(:,j));
                end
                obj.Y(:,i) = obj.solver(temp, dt, i);
            end
            
            % combine
            y = u0;
            for i = 1:obj.s
                y = y + dt*obj.b(i)*obj.L(dt + obj.c(i), obj.Y(:,i));
            end
            
            obj.u0 = y;
            obj.t = obj.t + dt;
        end
        
    end
    
    methods (Access = private)
        
        function  y = linearSolve(obj, y, dt,i)
            y = (obj.I - dt*obj.A(i,i)*obj.DT)\y;
        end
        
        function y = nonlinearImplicitStage( obj, y, dt, i )
            epss = 1e-14;
            options = optimset('Display','off', 'TolFun',epss, 'TolX',epss);
            fzero = @(K) (y + dt*obj.A(i,i)*obj.L(dt + obj.c(i), K))  - K;
            [y, ~] = fsolve(@(Y) fzero(Y), y, options);
        end
        
    end
    

end
