classdef BuckleyLeverett < handle
    
    properties
       name;
       f = @(t, u) u.^2./(u.^2 + 1.0/3.0*(1-u).^2);
       haveExact = false;
       eqn;
       CFL_MAX = 1;
       isLinear = false;
    end
    
    methods
        function obj = BuckleyLeverett(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.parse(varargin{:});
            
        end
    end
end