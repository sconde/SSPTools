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
        tFinal;
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
        BC;
    end
    
    methods
        function obj = LF(varargin)
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'isSSP', true);
            addParameter(p, 'Problem', []);
            addParameter(p, 'y0', []);
            addParameter(p, 't0', 0.0);
            addParameter(p, 'Tfinal', []);
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
            
            if ~isempty(p.Results.Tfinal)
                obj.tFinal = p.Results.Tfinal;
            end
            
            %TODO: is this true?
            obj.isSSP = true;
            
            obj.L = @(t,y) obj.problem.rhs(t, y);
            obj.dx = obj.problem.dx;
            obj.x = obj.problem.x(:); % make it a vector
            
            obj.systemSize = obj.problem.systemSize;
            obj.name = 'Lax?Friedrichs';
            
            % set the boundary condition function
            obj.BC = @(u, bcl, ul, bcr, ur) obj.problem.applyBC(u, bcl, ul, bcr, ur);
        end % LF constructor
    end
    
    methods %( Access = protected )
        
        function numFlux = numericalFlux(obj, maxvel, u, v)
            % function numFlux = numericalFlux(obj, u, v)
            
            % evaluate the flux
            fu = obj.problem.flux(u);
            fv = obj.problem.flux(v);
            numFlux = (fu + fv)/2 - maxvel/2*(v-u);
        end
        
        function FU = rhs(obj,maxvel, u) % acting like the L method?
            
            %Q = [r ru E];
            density   = u(1:obj.problem.N);
            momentum  = u(obj.problem.N+1:2*(obj.problem.N));
            energy    = u(2*(obj.problem.N)+1:3*(obj.problem.N));
            
            if strcmpi(obj.problem.ProblemType, 'sod')
                % Extend data and assign boundary conditions
                [~, re] = obj.BC(density, 'D', obj.problem.BCLvalue(1), 'D', obj.problem.BCRvalue(1));   % rho extended
                [~, me] = obj.BC(momentum, 'D', obj.problem.BCLvalue(2), 'N', obj.problem.BCRvalue(2));   % momentum extended
                [~, Ee] = obj.BC(energy, 'D', obj.problem.BCLvalue(3), 'N', obj.problem.BCRvalue(3));   % Energy extension
            elseif strcmpi(obj.problem.ProblemType, 'shock')
                [~, re] = obj.BC(density, 3.857143, 0);
                [~, me] = obj.BC(momentum, 10.141852, 0);
                [~, Ee] = obj.BC(energy, 39.1666661, 0);
            end
            
            Q = [re(2:obj.problem.N+1); me(2:obj.problem.N+1); Ee(2:obj.problem.N+1)];
            Qp = [re(3:obj.problem.N+2); me(3:obj.problem.N+2); Ee(3:obj.problem.N+2)];
            Qm = [re(1:obj.problem.N); me(1:obj.problem.N); Ee(1:obj.problem.N)];
            
            FU = -(obj.numericalFlux(maxvel, Q, Qp) - obj.numericalFlux(maxvel, Qm, Q))/obj.dx;
            
        end % rhs
        
        function [t] = takeStep(obj, dt)
            dt;
            % apply the boundary condition
            Q = obj.y0;

            [~, dt, maxvel] = obj.problem.closureModel(Q);
            
            if (obj.t0 + dt > obj.tFinal)
                dt = obj.tFinal - obj.t0;
            end
            dt;
            
            if ~isreal(dt)
                keyboard
            end
            
            lam = dt/obj.dx;
            FU = obj.rhs(maxvel, Q);
            maxvel
            
            u0 = Q + dt*FU;
            t = obj.t0 + dt;
            obj.y0 = u0;
            obj.t0 = t;
            
            if ~isreal(obj.t0)
                keyboard
            end
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
