classdef LinearAdvection < handle
    
    properties
       name;
       f;
       haveExact = false;
       eqn;
       y0;
       CFL_MAX;
       isLinear = true;
       a;
    end
    
    methods
        function obj = LinearAdvection(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            addParamValue(p, 'y0', []);
            addParamValue(p, 'a', 1);
            p.parse(varargin{:});
            if isa(p.Results.y0, 'function_handle')
                obj.y0 = p.Results.y0;
                obj.haveExact = true;
            else
                obj.y0 = p.Results.y0(:);
            end
            
            obj.a = p.Results.a;
            obj.CFL_MAX = abs(obj.a);
            obj.f = @(t, u) obj.a *u;
        end
    end
end