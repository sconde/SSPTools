classdef ERK < SSPTools.Steppers.RK
    
    properties
        
    end
    
    properties ( Access = protected)
        isExplicit = true;
        isMSRK = true; 	% Multi-Stage Runge-Kutta
        isButcher = true; % starting with Butcher formulation
        isLowStorage = false;  % need a way to determine is low-storage
        n;
        Y;
        Fvec;
    end
    
    methods
        
        function obj = ERK(varargin)
            obj = obj@SSPTools.Steppers.RK(varargin{:});
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'name','MSRK-ERK');
            addParameter(p, 'isSSP', false);
            addParameter(p, 'isButcher', true);
            addParameter(p, 'isLowStorage', false);
            p.parse(varargin{:});

            if isa(obj.y0, 'function_handle')
                obj.u0 = obj.y0(obj.x);
                obj.n = size(obj.x,1);
            else
                obj.u0 = obj.y0;
                obj.n = size(obj.u0,1);
            end
            
            obj.name = p.Results.name;
            obj.Y = zeros(obj.n, obj.s);
            obj.Fvec = zeros(obj.n, obj.s);
            
            if isa(obj.dfdx, 'WenoCore.Weno')
                obj.dfdx.f = @(t,u) obj.ExplicitProblem.f(t,u);
                obj.dfdx.em = obj.ExplicitProblem.em;
            end
            
            obj.verifyMethod();
        end % end constructor
        
        
    end
    
    methods %( Access = protected )
        
        function [y, dt] = takeStep(obj, dt)
            % function [y, dt] = takeStep(dt)
            % returns the new solution (y) and time-step taken (dt)
            u0 = obj.u0;
            t_ = obj.t;
            obj.Y(:,1) = u0;
            obj.Fvec(:,1) = obj.L(t_ , u0);
            obj.dt_ = dt;
            
            % intermediate stage value
            for i = 2:obj.s
                temp = u0;
                for j = 1:i-1
                    temp = temp + dt*obj.A(i,j)*obj.Fvec(:, j);
                end
                obj.Fvec(:, i) = obj.L(t_ + obj.c(i)*dt, temp);
            end
            
            % combine
            y = u0 + dt*obj.Fvec*obj.b(:);
            obj.u0 = y;
            obj.t  = obj.t + dt;
        end
        
    end
    
    methods (Access = private)
        
        function verifyMethod(obj)
            assert(isequal(tril(obj.A,-1),obj.A),...
                'Method is not Explicit');
        end
        
    end
    

end
