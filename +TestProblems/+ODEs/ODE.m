classdef ODE < handle
    
    properties
        name;
        f;
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
