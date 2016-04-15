classdef PDE < handle
    
    properties
       name;
       f; %flux
       haveExact = false;
       eqn;
       y0;
       CFL_MAX = 1;
    end
    
    methods
        function obj = PDE(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            addParamValue(p, 'f', []);
            addParamValue(p, 'y0', []);
            p.parse(varargin{:});
            
            if ~isempty(p.Results.f)
                obj.f = p.Results.f;
            end
            
            if isa(p.Results.y0, 'function_handle')
                obj.y0 = p.Results.y0;
                obj.haveExact = true;
            else
                obj.y0 = p.Results.y0(:);
            end
        end
    end
end