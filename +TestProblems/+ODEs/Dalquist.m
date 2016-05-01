classdef Dalquist < TestProblems.ODEs.ODE
    
    properties
        a;
        isLinear;
        exact;
    end
    
    methods
        
        function obj = Dalquist(varargin)
            
            obj = obj@TestProblems.ODEs.ODE(varargin{:});
            p = inputParser;
            p.KeepUnmatched = true;
            
            
            addParameter(p, 'a', 2);
            p.parse(varargin{:});
            
            obj.name = 'Dalquist';
            obj.isLinear = false;
            
            obj.a = p.Results.a;
            obj.f = @(t, u) obj.a*u;
            obj.exact = @(t,u) exp(obj.a*t);
        end
        
    end
    
    methods (Access = protected)
        function dudt = rhs(ep, u)
            dudt = [u(2);(1/ep)*(-u(1) + (1 - u(1)^2)*u(2))];
        end
    end
end
