classdef AdditiveODE < handle
    
    properties
        name;
        f;
        g;
        isSystem = false;
    end
    
    methods
        
        function obj = ODE(varargin)
            
        end
        
    end
    
    methods (Access = protected)
        function dudt = Frhs(u) % explicit function
        end
        
        function dudt = Grhs(u) % implicit function
        end
    end
end
