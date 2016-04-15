classdef LinearAdvection < TestProblems.PDES.PDE
    
    properties
        a;
    end
    
    methods
        function obj = LinearAdvection(varargin)
            obj = obj@TestProblems.PDES.PDE(varargin{:});
            
            keyboard
            p = inputParser;
            p.KeepUnmatched = true;
            addParamValue(p, 'a', 1);
            addParamValue(p, 'y0', []);
            p.parse(varargin{:});
            
            obj.a = p.Results.a;
            obj.y0 = p.Results.y0;
            obj.flux = @(u) obj.a * u;
            obj.name = 'Linear-Advection';
            obj.eqn = 'u_t + a*(u)_x = 0';
            
            if isa(obj.y0, 'function_handle')
                obj.haveExact = true;
            else
                obj.y0 = obj.y0(:);
            end
        end
        
    end
end