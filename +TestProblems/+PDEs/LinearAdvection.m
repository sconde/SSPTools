classdef LinearAdvection < handle
    
    properties
       name = 'Linear Advection';
       f;
       haveExact = false;
       eqn;
       CFL_MAX;
       isLinear = true;
       a;
    end
    
    methods
        function obj = LinearAdvection(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            
            addParamValue(p, 'a', 1);
            p.parse(varargin{:});

            
            obj.a = p.Results.a;
            obj.CFL_MAX = abs(obj.a);
            obj.f = @(t, u) obj.a *u;
        end
    end
end