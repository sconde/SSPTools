classdef Vanderpol < TestProblems.ODEs.ODE
    
    properties
        ep;
        isLinear;
    end
    
    methods
        
        function obj = Vanderpol(varargin)
            
            obj = obj@TestProblems.ODEs.ODE(varargin{:});
            p = inputParser;
            p.KeepUnmatched = true;
            
            
            addParameter(p, 'epsilon', 0.1);
            p.parse(varargin{:});
            
            obj.name = 'Vanderpol';
            obj.isLinear = false;
            
            obj.ep = p.Results.epsilon;
            %obj.f = @(t, y) rhs(obj.ep, y);
            obj.f = @(t, u) [u(2);(1/obj.ep)*(-u(1) + (1 - u(1)^2)*u(2))];
        end
        
    end
    
    methods (Access = protected)
        function dudt = rhs(ep, u)
            dudt = [u(2);(1/ep)*(-u(1) + (1 - u(1)^2)*u(2))];
        end
    end
end
