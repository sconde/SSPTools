classdef LinearAdvection < handle
    
    properties
       name = 'Linear Advection';
       f;
       haveExact = false;
       eqn;
       CFL_MAX;
       isLinear = true;
       a;
       em;
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
            obj.em = @(u) obj.a*ones(size(u));
        end
    end
end