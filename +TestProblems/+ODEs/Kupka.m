classdef Kupka < TestProblems.ODEs.AdditiveODE
    
    properties
        a;
        isLinear;
        exact;
        isSystem = false;
    end
    
    methods
        
        function obj = Kupka(varargin)
            
            obj = obj@TestProblems.ODEs.ODE(varargin{:});
            p = inputParser;
            p.KeepUnmatched = true;
            
            
            addParameter(p, 'a', 2);
            p.parse(varargin{:});
            
            obj.name = 'Dalquist';
            obj.isLinear = false;
            
            obj.a = p.Results.a;
            obj.f = @(t, u) ( 1 + sin(u));
            obj.g = @(t, u) (u.^2 - sin(u));
            obj.exact = @(t,u) tan(t);
        end
        
    end
    
    methods (Access = protected)
        function dudt = Frhs(ep, u)
            dudt = [u(2);(1/ep)*(-u(1) + (1 - u(1)^2)*u(2))];
        end
    end
end
