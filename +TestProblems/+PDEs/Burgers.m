classdef Burgers < handle
    
    properties
       name = 'Burgers';
       f = @(t, u) 0.5*(u.^2);
       haveExact = false;
       eqn;
       CFL_MAX = 1;
       isLinear = false;
       em;
       isSystem = false;
       systemSize = 1;
       xx;
    end
    
    methods
        function obj = Burgers(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            addParamValue(p, 'y0', []);
            p.parse(varargin{:});
            
            obj.em = @(u) u;
            
        end
        
        function set.xx(obj, fin)
            obj.xx = fin;
        end
    end
end
