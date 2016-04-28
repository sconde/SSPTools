classdef Weno < SSPTools.Discretizers.Discretize
    
    properties
        isSSP = true;
        name = 'Finite Difference';
        bc;
        nx;
    end
    
    methods
        
        function obj = Weno(varargin)
            obj = obj@SSPTools.Discretizers.Discretize();
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'x', []);
            addParameter(p, 'isSSP', true);
            addParameter(p, 'derivativeOrder', 1);
            addParameter(p, 'bc', 'periodic');
            p.parse(varargin{:});
            
            obj.x = p.Results.x;
            obj.isSSP = p.Results.isSSP;
            obj.bc = p.Results.bc;
            obj.derivativeOrder = p.Results.derivativeOrder;
            
            if ~isempty(obj.x)
                obj.nx = numel(obj.x);
                obj.x = obj.x(:);
            end
        end
        
    end
    
    methods
        function [y_x] = L(obj, y) end
    end
    
    
end
