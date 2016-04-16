classdef Burgers < handle
    
    properties
       name;
       f = @(t, u) 0.5*(u.^2);
       haveExact = false;
       eqn;
       y0;
       CFL_MAX = 1;
       isLinear = false;
    end
    
    methods
        function obj = Burgers(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            addParamValue(p, 'y0', []);
            p.parse(varargin{:});
            
            if isa(p.Results.y0, 'function_handle')
                obj.y0 = p.Results.y0;
                obj.haveExact = true;
            else
                obj.y0 = p.Results.y0(:);
            end
        end
    end
end