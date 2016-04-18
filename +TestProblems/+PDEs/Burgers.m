classdef Burgers < handle
    
    properties
       name;
       f = @(t, u) 0.5*(u.^2);
       haveExact = false;
       eqn;
       CFL_MAX = 1;
       isLinear = false;
    end
    
    methods
        function obj = Burgers(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            addParamValue(p, 'y0', []);
            p.parse(varargin{:});
            
        end
    end
end