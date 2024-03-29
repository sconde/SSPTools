classdef CurtissHirschfelder < TestProblems.ODEs.ODE
% A First Course in the Numerical Analysis of Differential Equations (Second Edition)
% Arieh Iserles, University of Cambridge
% Page 107, equation 6.4

    properties
        isLinear;
        exact;
    end
    
    methods
        
        function obj = CurtissHirschfelder(varargin)
            
            obj = obj@TestProblems.ODEs.ODE(varargin{:});
            p = inputParser;
            p.KeepUnmatched = true;
            
            
            obj.name = 'Curtiss-Hirschfelder';
            obj.isLinear = false;
            
            obj.f = @(t, u) -50*(u - cos(t))
            obj.exact = @(t,u) (2500/2501)*cos(t) + (50/2501)*sin(t) + (1/2501)*exp(-50*t);
        end
        
    end
    
end
