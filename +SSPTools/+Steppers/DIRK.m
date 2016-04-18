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
        u0;
    end
    
    methods
        
        function obj = DIRK(varargin)
            obj = obj@SSPTools.Steppers.RK(varargin{:});
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParamValue(p,'name','MSRK-DIRK');
            addParamValue(p, 'isSSP', false);
            addParamValue(p, 'isButcher', true);
            addParamValue(p, 'isLowStorage', false);
            p.parse(varargin{:});

            obj.name = p.Results.name;  
            obj.n = size(obj.x,1);
            obj.Y = zeros(obj.n, obj.s);
            obj.I = speye(obj.n);
            obj.isImplicitLinear = obj.ExplicitProblem.isLinear;

            if obj.isImplicitLinear
                obj.solver = @(y, dt, i) linearSolve(obj, y, dt,i);
                obj.DT = obj.dfdx.D;
                obj.NL = @(t,y) obj.dfdx.L(obj.ExplicitProblem.f(t,y));
            else
                obj.NL = @nonlinearImplicitStage;
            end
            
            if isa(obj.y0, 'function_handle')
                obj.u0 = obj.y0(obj.x);
            else
                obj.u0 = obj.y0;
            end
        end
        
        
    end
    
    methods %( Access = protected )
        function [y] = takeStep(obj, dt)
            
            %check to see if CFL violation
            assert((dt/obj.dx) <= obj.CFL, ...
                sprintf('ERK: CFL Violation (CFL = %3.2f )',dt/obj.dx) );
            
            u0 = obj.u0;
            
            % first stage implicit solve
            obj.Y(:,1) = obj.solver(u0,dt, 1);
            
            % intermediate stage value
            for i = 1:obj.s-1
                
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
            
        end
        
    end
    
    methods (Access = private)
        
        function  y = linearSolve(obj, y, dt,i)
            y = (obj.I - dt*obj.A(i,i)*obj.DT)\y;
        end
        
        function y = nonlinearImplicitStage( y )
            
        end
        
    end
    

end
