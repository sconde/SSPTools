classdef RK < handle
    
    properties
        name;			% Name of time-stepping method.
        dfdt;
        dfdx;
        ExplicitProblem;
        A = []; b = []; c = []; alpha = []; s = []; r = [];
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
            addParamValue(p,'name','MSRK');
            addParamValue(p, 'dfdx', []);
            addParamValue(p, 'dfdt', []);
            addParamValue(p, 'ExplicitProblem', []);
            addParamValue(p, 'ImplicitProblem', []);
            addParamValue(p, 'A', []);
            addParamValue(p, 'b', []);
            addParamValue(p, 's', []);
            addParamValue(p, 'r', []);
            addParamValue(p, 'alpha', []);
            addParamValue(p, 't', 0.0);
            addParamValue(p, 't0', 0);
            addParamValue(p, 'y0', []);
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
            obj.s = p.Results.s;
            obj.alpha = p.Results.alpha;
            obj.name = p.Results.name;
            obj.dfdx = p.Results.dfdx;
            obj.t0 = p.Results.t0;
            assert(isequal(p.Results.s, size(p.Results.A,1)),...
                sprintf('RK A:Stage-count -- Num-Rows(A) != %d',p.Results.s));
            
            if ~isempty(p.Results.r)
                obj.r = p.Results.r;
                obj.isSSP = true;
            end
            
            if ~isempty(p.Results.ExplicitProblem)
                obj.ExplicitProblem = p.Results.ExplicitProblem;
            end
            
            obj.L = @(t,y) obj.dfdx.L(obj.ExplicitProblem.f(t, y));
            
%             if isa(obj.ExplicitProblem.y0, 'function_handle')
%                 obj.y0 = obj.ExplicitProblem.y0(obj.dfdx.x);
%             else
%                 obj.y0 = obj.ExplicitProblem.y0(:);
%             end
            
            obj.CFL = obj.ExplicitProblem.CFL_MAX;
            obj.dx = obj.dfdx.dx;
            obj.x = p.Results.dfdx.x;
            obj.isLinear = obj.ExplicitProblem.isLinear;
        end
    end
    
    methods %( Access = protected )
        function [y] = takeStep(obj, dt) end
        
        function resetInitCondition(obj)
            obj.u0 = obj.y0(obj.x);
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
