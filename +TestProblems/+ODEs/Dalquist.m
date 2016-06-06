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
            obj.isLinear = true;
            
            obj.a = p.Results.a;
            obj.f = @(t, u) obj.a*u;
            obj.exact = @(t,u) exp(obj.a*t);
        end
        
    end
    
end
