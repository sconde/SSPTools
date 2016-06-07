classdef Brusselator < TestProblems.ODEs.ODE
    
    properties
        ep;
        isLinear;
    end
    
    methods
        
        function obj = Brusselator(varargin)
            
            obj = obj@TestProblems.ODEs.ODE(varargin{:});
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.parse(varargin{:});
            
            obj.name = 'Brusselator';
            obj.isLinear = false;
            
            obj.f = @(t, u) [ 1 + u(2)*u(1)^2 - 4*u(1);
                3*u(1) - u(2)*u(1)^2];
        end
        
    end
    
    methods (Access = protected)

    end
end
