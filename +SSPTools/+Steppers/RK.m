classdef RK < handle
    
    properties
        name;			% Name of time-stepping method.
        dfdt;
        dfdx;
        ExplicitProblem;
        A = []; b = []; c = []; alpha = []; s = []; r = [];
        p; % Nonlinear-Order
        plin; % Linear-Order
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
        u0;
    end
    
    methods
        function obj = RK(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'name','MSRK');
            addParameter(p, 'dfdx', []);
            addParameter(p, 'dfdt', []);
%             addParameter(p, 'ExplicitProblem', []);
%             addParameter(p, 'ImplicitProblem', []);
            addParameter(p, 'A', []);
            addParameter(p, 'b', []);
            addParameter(p, 'r', []);
            addParameter(p, 'p', []);
            addParameter(p, 'plin', []);
            addParameter(p, 'alpha', []);
            addParameter(p, 't', 0.0);
            addParameter(p, 't0', 0);
            addParameter(p, 'y0', []);
            p.parse(varargin{:});
            
            if isa(p.Results.dfdt, 'function_handle')
                obj.dydt = p.Results.dfdt;
            end
            
            if isa(p.Results.y0, 'function_handle')
                obj.y0 = p.Results.y0;
                obj.haveExact = true;
            else
                obj.y0 = p.Results.y0(:);
            end
                        
            obj.t = p.Results.t;
            obj.A = p.Results.A;
            obj.b = p.Results.b;
            obj.c = sum(obj.A,2);
            obj.s = numel(obj.b); %infer the number of stages from the length of b
            obj.p = p.Results.p;
            
            if isempty(p.Results.plin)
                obj.plin = obj.p;
            else
                obj.plin = p.Results.plin;
            end

            obj.alpha = p.Results.alpha;
            obj.name = p.Results.name;
            obj.dfdx = p.Results.dfdx;
            obj.t0 = p.Results.t0;
            
            assert(isequal(obj.s, size(p.Results.A,1)),...
                sprintf('RK A:Stage-count -- Num-Rows(A) != %d',obj.s));
            
            if ~isempty(p.Results.r)
                obj.r = p.Results.r;
                obj.isSSP = true;
            end
                        
            if ~isempty(obj.dfdx.problem)
                obj.ExplicitProblem = obj.dfdx.problem;
            end
                        
            if isa(obj.ExplicitProblem, 'TestProblems.ODEs.ODE')
                obj.L = @(t,y) obj.ExplicitProblem.f(t,y);
                obj.CFL = [];
                obj.dx = [];
                obj.x = [];
            else
                if ~isa(obj.dfdx, 'WenoCore.Weno')
                    obj.L = @(t,y) obj.dfdx.L(obj.ExplicitProblem.f(t, y));
                else
                    obj.L = @(t,y) obj.dfdx.L(t, y);
                end
                
                %obj.CFL = obj.ExplicitProblem.CFL_MAX;
                obj.dx = obj.dfdx.dx;
                obj.x = p.Results.dfdx.x;
            end

            obj.isLinear = obj.ExplicitProblem.isLinear;
            
            if obj.isSSP
                obj.name = sprintf('SSP(%d,%d)%d',obj.s, obj.p, obj.plin);
            else
                obj.name = sprintf('RK(%d,%d)%d',obj.s, obj.p, obj.plin);
            end
            
        end
    end
    
    methods %( Access = protected )
        function [y] = takeStep(obj, dt) end
        
        function resetInitCondition(obj)
            if isa(obj.y0, 'function_handle')
                obj.u0 = obj.y0(obj.x);
            else
                obj.u0 = obj.y0;
            end
            obj.t = 0.0;
        end
        
        function [t, y] = getState(obj)
            y = obj.u0;
            t = obj.t;
        end
        
    end
    
    methods ( Access = private )
        
        function obj = setL(obj)
            obj.L = @(t, y) obj.L(y);
        end
        
    end
    
end
