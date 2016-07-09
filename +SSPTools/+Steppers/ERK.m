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

            
            if obj.ExplicitProblem.isSystem
                U1 = obj.y0{1}(obj.dfdx.x, obj.dfdx.y, 0);
                U2 = obj.y0{2}(obj.dfdx.x, obj.dfdx.y, 0);
                obj.u0 = [U1(:); U2(:)];
                obj.n = size(obj.u0,1);
                
                assert(isequal(obj.ExplicitProblem.systemSize, 2),...
                    'CodeLimitation: Only working for 2D at the moment');
            else
                
                if isa(obj.y0, 'function_handle')
                    obj.u0 = obj.y0(obj.x);
                    obj.n = size(obj.x,1);
                else
                    obj.u0 = obj.y0;
                    obj.n = size(obj.u0,1);
                end
            end
            
            obj.name = p.Results.name;
            obj.Y = zeros(obj.n, obj.s);
            obj.Fvec = zeros(obj.n, obj.s);
            
            if isa(obj.dfdx, 'WenoCore.Weno')
                                
                obj.dfdx.f = obj.ExplicitProblem.f;
                obj.dfdx.em = obj.ExplicitProblem.em;
            elseif ( obj.ExplicitProblem.systemSize > 1)
                obj.L = @(t,u) obj.ExplicitProblem.f(t,u);
            end
            
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
    

end
