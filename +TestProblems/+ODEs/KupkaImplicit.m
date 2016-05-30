classdef KupkaImplicit < TestProblems.ODEs.ODE
    
    properties
        a;
        isLinear;
        exact;
    end
    
    methods
        
        function obj = KupkaImplicit(varargin)
            
            obj = obj@TestProblems.ODEs.ODE(varargin{:});
            p = inputParser;
            p.KeepUnmatched = true;
            
            
            addParameter(p, 'a', 2);
            p.parse(varargin{:});
            
            obj.name = 'Kupka-Explicit';
            obj.isLinear = false;
            
            obj.a = p.Results.a;
            obj.f = @(t, u) (u.^2 - sin(u));
            obj.exact = @(t,u) zero(size(u));
        end
        
    end
    
    methods (Access = protected)
        function dudt = rhs(ep, u)
            dudt = [u(2);(1/ep)*(-u(1) + (1 - u(1)^2)*u(2))];
        end
    end
end
