classdef ODE < handle
    
    properties
        name;
        f;
        isSystem = false;
        systemSize = 1;
    end
    
    methods
        
        function obj = ODE(varargin)
            
        end
        
    end
    
    methods (Access = protected)
        function dudt = rhs(u)
        end
    end
end
