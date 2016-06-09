classdef ERK < SSPTools.Steppers.RK
    
    properties
        
    end
    
    properties ( Access = private)
        isExplicit = true;
        isMSRK = true; 	% Multi-Stage Runge-Kutta
        isButcher = true; % starting with Butcher formulation
        isLowStorage = false;  % need a way to determine is low-storage
        n;
        Y;
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
            
            if isa(obj.dfdx, 'WenoCore.Weno')
                obj.dfdx.f = obj.ExplicitProblem.f;
                obj.dfdx.em = obj.ExplicitProblem.em;
            end
        end % end constructor
        
        
    end
    
    methods %( Access = protected )
        
        function [y, dt] = takeStep(obj, dt)
            % function [y, dt] = takeStep(dt)
            % returns the new solution (y) and time-step taken (dt)
            u0 = obj.u0;
            obj.Y(:,1) = u0;
            obj.dt_ = dt;
            
            % intermediate stage value
            for i = 2:obj.s
                temp = u0;
                for j = 1:i-1
                    %keyboard
                    temp = temp + dt*obj.A(i,j)*obj.L(dt + obj.c(j), obj.Y(:,j));
                    %keyboard
                end
                obj.Y(:,i) = temp;
            end
            
            % combine
            y = u0;
            for i = 1:obj.s
                y = y + dt*obj.b(i)*obj.L(dt + obj.c(i), obj.Y(:,i));
            end
            
            obj.u0 = y;
            obj.t  = obj.t + dt;
        end
        
    end
    

end
