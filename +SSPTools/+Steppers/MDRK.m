classdef MDRK < handle
    % MDRK stepper
    
    properties
        name;			% Name of time-stepping method.
        y0;
    end
    
    properties (Access = protected)
        steps = 1;
        isSSP;
        L;
        CFL;
        dx;
        isLinear;
    end
    
    methods
        function obj = MDRK(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParamValue(p,'name','BDf');
            addParamValue(p, 't0', 0);
            p.parse(varargin{:});
            
            if isa(p.Results.dfdt, 'function_handle')
                obj.dydt = p.Results.dfdt;
            end
            
            
            obj.name = p.Results.name;
            obj.dfdx = p.Results.dfdx;
            obj.t0 = p.Results.t0;
            
        end
    end
    
    methods
        function [y] = takeStep(obj, dt) end
        
    end
    
    methods ( Access = private )
        
        function obj = setL(obj)
            obj.L = @(t, y) obj.L(y);
        end
        
    end
    
end