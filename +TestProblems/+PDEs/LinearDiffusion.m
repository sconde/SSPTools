classdef LinearDiffusion < handle
    
    properties
       name = 'Linear-Diffusion';
       f;
       haveExact = false;
       eqn;
       isLinear = true;
       nu;
       isSystem = false;
       systemSize = 1;
    end
    
    methods
        function obj = LinearDiffusion(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            addParameter(p, 'nu', 0.01);
            p.parse(varargin{:});

            
            obj.nu = p.Results.nu;
            obj.f = @(t, u) -obj.nu *u;
        end
    end
end
