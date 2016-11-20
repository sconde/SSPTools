classdef LinearAdvection < handle
    
    properties
       name = 'Linear-Advection';
       f;
       haveExact = false;
       eqn;
       CFL_MAX;
       isLinear = true;
       a;
       em;
       isSystem = false;
       systemSize = 1;
       xx;
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

        function set.xx(obj, fin)
            obj.xx = fin;
        end
    end
end
