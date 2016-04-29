classdef IMEXRK < SSPTools.Steppers.RK
    
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
        utex;
        utimp;
    end
    
    methods
        
        function obj = IMEXRK(varargin)
            obj = obj@SSPTools.Steppers.RK(varargin{:});
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'name','MSRK-IMEXRK');
            addParameter(p, 'isSSP', false);
            addParameter(p, 'isButcher', true);
            addParameter(p, 'ImplicitProblem', []);
            addParameter(p, 'At', []);
            addParameter(p, 'bt', []);
            addParameter(p, 'isLowStorage', false);
            addParameter(p, 'dgdx', []);
            addParameter(p, 't', 0.0);
            p.parse(varargin{:});
            
            obj.At = p.Results.At;
            obj.bt = p.Results.bt;
            obj.ct = sum(obj.A,2);
            obj.name = p.Results.name;
            obj.n = size(obj.x,1);
            obj.I = speye(obj.n);
            obj.Y = zeros(obj.n, obj.s);
            obj.F = zeros(obj.n, obj.s);
            obj.G = zeros(obj.n, obj.s);
            obj.t = p.Results.t;
            
            if ~isempty(p.Results.ImplicitProblem)
                obj.ImplicitProblem = p.Results.ImplicitProblem;
                obj.isImplicitLinear = p.Results.ImplicitProblem.isLinear;
            end
            
            if ~isempty(p.Results.dgdx)
                obj.dgdx = p.Results.dgdx;
            end
            
            if obj.isImplicitLinear
                obj.solver = @(y, dt, i) linearSolve(obj, y, dt,i);
                obj.F = @(t,y) obj.dfdx.L(obj.ExplicitProblem.f(t,y));
                obj.G = @(t,y) obj.dgdx.L(obj.ImplicitProblem.f(t,y));
                obj.DT = obj.dgdx.D;
            else
                obj.solver = @(y, dt, i) nonlinearImplicitStage( obj, y, dt, i );
                obj.F = @(t,y) obj.dfdx.L(obj.ExplicitProblem.f(t,y));
                obj.G = @(t,y) obj.dgdx.L(obj.ImplicitProblem.f(t,y));
            end
            
            if isa(obj.y0, 'function_handle')
                obj.u0 = obj.y0(obj.x);
            else
                obj.u0 = obj.y0;
            end
        end
        
        
    end
    
    methods


        function [y] = takeStep(obj, dt)
            
%             %check to see if CFL violation
%             assert((dt/obj.dx) <= obj.CFL, ...
%                 sprintf('ERK: CFL Violation (CFL = %3.2f )',dt/obj.dx) );
%             
            u0 = obj.u0;
            
            % first stage implicit solve
            obj.Y(:,1) = obj.solver(u0,dt, 1);
            
            % intermediate stage value
            for i = 2:obj.s
                
                obj.utex(:,i) = obj.F(dt + obj.c(i), obj.Y(:,i));
                obj.utimp(:,i) = obj.G(dt + obj.ct(i), obj.Y(:,i));
                
                temp = u0;
                for j = 1:i-1
                    temp = temp + dt*obj.A(i,j)*obj.F(dt + obj.c(j), obj.Y(:,j)) + ...
                        dt*obj.At(i,j)*obj.G(dt + obj.ct(j), obj.Y(:,j));
                end
                obj.Y(:,i) = obj.solver(temp, dt, i);
            end
            
            % combine
            y = u0;
            for i = 1:obj.s
                y = y + dt*obj.b(i)*obj.F(dt + obj.c(i), obj.Y(:,i)) +...
                    dt*obj.bt(i)*obj.G(dt + obj.c(i), obj.Y(:,i));
            end
            
            obj.u0 = y;
            obj.t = obj.t + dt;
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
