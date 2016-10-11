classdef Mathiew < TestProblems.ODEs.ODE
    
    properties
        a;b;
        isLinear;
    end
    
    methods
        
        function obj = Mathiew(varargin)
            
            obj = obj@TestProblems.ODEs.ODE(varargin{:});
            p = inputParser;
            p.KeepUnmatched = true;
            
            
            addParameter(p, 'a', 2);
            addParameter(p, 'b', 1);
            p.parse(varargin{:});
            
            obj.name = 'Mathiew';
            obj.isLinear = false;
            
            obj.a = p.Results.a;
            obj.b = p.Results.b;
            obj.f = @(t, u) [u(2); -(obj.a - obj.b*cos(2*t))*u(1)];
        end
        
    end
    
    methods (Access = protected)

    end
end
