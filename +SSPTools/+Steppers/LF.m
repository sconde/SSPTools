classdef LF < handle
    
    properties
        name;			% Name of time-stepping method.
        dfdt;
        dfdx;
        problem;
        p; % Order
        t0;
        y0;
        t;
    end
    
    properties (Access = protected)
        steps = 1;
        isSSP;
        L;
        CFL;
        dx;
        isLinear;
        haveExact = false;
        x;
        systemSize;
    end
    
    methods
        function obj = LF(varargin)
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'isSSP', true);
            addParameter(p, 'Problem', []);
            addParameter(p, 'y0', []);
            addParameter(p, 't0', 0.0);
            p.parse(varargin{:});
            
            p.parse(varargin{:});
            
            % get the problem definition (i.e flux)
            if ~isempty(p.Results.Problem)
                obj.problem = p.Results.Problem;
            elseif ~isempty(p.Results.ODE)
                obj.problem = p.Results.ODE;
            end
            
            if isa(p.Results.y0, 'function_handle')
                obj.y0 = p.Results.y0;
                obj.haveExact = true;
            else
                obj.y0 = p.Results.y0(:);
            end
            
            obj.t0 = p.Results.t0;
            
            %TODO: is this true?
            obj.isSSP = true;
            
            obj.L = @(t,y) obj.problem.rhs(t, y);
            obj.dx = obj.problem.dx;
            obj.x = obj.problem.x(:); % make it a vector
            
            obj.systemSize = obj.problem.systemSize;
            obj.name = 'Lax?Friedrichs';
            
        end % LF constructor
    end
    
    methods %( Access = protected )
        
        function numFlux = numericalFlux(obj, dt, t, u, v)
           % function numFlux = numericalFlux(obj, u, v)
           
           lam = dt/obj.dx;
           % evaluate the flux
           fu = obj.ExplicitProblem.f(t, u);
           fv = obj.ExplicitProblem.f(t, v);
           numFlux = (fu + fv)/2 - lam/2*(v-u);
        end
        
        function [t0, u0] = takeStep(obj, dt)
            Q = obj.y0;
            
            % get the time scale
            [~, dt] = obj.problem.closureModel(Q);
 
            u0 = Q + dt*obj.L(dt, Q);
            obj.y0 = u0;
            obj.t0 = obj.t0 + dt;
        end
        
        function resetInitCondition(obj)
            if isa(obj.y0, 'function_handle')
                obj.u0 = obj.y0(obj.x);
            else
                obj.u0 = obj.y0;
            end
            obj.t = 0.0;
        end
        
        function [t, y] = getState(obj)
            y = obj.y0;
            t = obj.t0;
            
            if (obj.systemSize > 1)
                y = reshape(y, obj.problem.N, obj.systemSize);
            end
            
        end
        
    end
    
    methods ( Access = private )
        
        function obj = setL(obj)
            obj.L = @(t, y) obj.L(y);
        end
        
    end
    
end
