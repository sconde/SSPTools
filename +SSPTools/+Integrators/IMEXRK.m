classdef IMEXRK < SSPTools.Integrators.RK
    
    properties
        At; bt; ct; dgdx; dgdt; ImplicitProblem;

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
        isImplicitLinear = true; %so far handling linear advection
        NL;
        DT;
        solver;
        I;
    end
    
    methods
        
        function obj = IMEXRK(varargin)
            obj = obj@SSPTools.Integrators.RK(varargin{:});
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParamValue(p,'name','MSRK-IMEXRK');
            addParamValue(p, 'isSSP', false);
            addParamValue(p, 'isButcher', true);
            addParamValue(p, 'ImplicitProblem', []);
            addParamValue(p, 'At', []);
            addParamValue(p, 'bt', []);
            addParamValue(p, 'isLowStorage', false);
            addParamValue(p, 'dgdx', []);
            p.parse(varargin{:});
            
            obj.At = p.Results.At;
            obj.bt = p.Results.bt;
            obj.ct = sum(obj.A,2);
            obj.name = p.Results.name;
            obj.n = size(obj.y0,1);
            obj.I = speye(obj.n);
            obj.Y = zeros(obj.n, obj.s);
            obj.F = zeros(obj.n, obj.s);
            obj.G = zeros(obj.n, obj.s);
            
            if ~isempty(p.Results.ImplicitProblem)
                obj.ImplicitProblem = p.Results.ImplicitProblem;
                obj.isImplicitLinear = p.Results.ImplicitProblem.isLinear;
            end
            
            if ~isempty(p.Results.dgdx)
                obj.dgdx = p.Results.dgdx;
            end
                        
            if obj.isImplicitLinear
                obj.solver = @(y, dt, i) linearSolve(obj, y, dt,i);
                obj.DT = obj.dgdx.D;
                obj.NL = @(t,y) obj.dgdx.L(obj.ImplicitProblem.f(t,y));
            else
                obj.NL = @nonlinearImplicitStage;
            end
        end
        
        
    end
    
    methods
        function [y] = takeStep(obj, dt)
            
            %check to see if CFL violation
            assert((dt/obj.dx) <= obj.CFL, ...
                sprintf('ERK: CFL Violation (CFL = %3.2f )',dt/obj.dx) );
            
            u0 = obj.y0;
            
            % first stage implicit solve
            obj.Y(:,1) = obj.solver(u0,dt, 1);
            
            
            % intermediate stage value
            for i = 1:obj.s
                
                obj.F(:,i) = obj.dfdx.L((obj.ExplicitProblem.f(dt + obj.c(i), obj.Y(:,i))));
                obj.G(:,i) = obj.dgdx.L(obj.ImplicitProblem.f(dt + obj.ct(i), obj.Y(:,i)));
                
                temp = u0;
                for j = 1:i-1
                    temp = temp + dt*obj.A(i,i)*obj.L(dt + obj.c(j), obj.Y(:,j)) + ...
                        dt*obj.At(i,j)*obj.NL(dt + obj.ct(j), obj.Y(:,j));
                end
                obj.Y(:,i) = obj.solver(u0,dt, i);
            end
            
            % combine
            y = u0;
            for i = 1:obj.s
                y = y + dt*obj.b(i)*obj.L(dt + obj.c(i), obj.Y(:,i)) +...
                    dt*obj.bt(i)*obj.NL(dt + obj.ct(i), obj.Y(:,i));
            end
            
            obj.y0 = y;
            
        end
        
    end
    
    methods (Access = private)
        
        function  y = linearSolve(obj, y, dt,i)
            y = (obj.I - dt*obj.At(i,i)*obj.DT)\y;
        end
        
        function y = nonlinearImplicitStage( y )
            
        end
        
    end
    
    
end
