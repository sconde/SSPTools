classdef Robertson < TestProblems.ODEs.ODE
    
    properties
        isLinear;
    end
    
    methods
        
        function obj = Robertson(varargin)
            
            obj = obj@TestProblems.ODEs.ODE(varargin{:});
            p = inputParser;
            p.KeepUnmatched = true;
            
            
            p.parse(varargin{:});
            
            obj.name = 'Chemical Reaction of Robertson';
            obj.isLinear = false;
            
            obj.f = @(t, u) [-0.04*u(1) + (10^4)* u(2)*u(3);
                0.04*u(1) - (10^4)*u(1)*u(3) - 3*10^7 *u(2);
                3*10^7*u(2)];
        end
        
    end
    
    methods (Access = protected)

    end
end
