classdef LoadERK < SSPTools.Steppers.ERK
    
    properties
        MethodName;
    end
    
    properties (Access = private)
        gam;
    end
    
    methods
        function obj = LoadERK(varargin)
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'MethodName','FE');
            addParameter(p, 'Theta', 0.1);
            p.parse(varargin{:});
            
            gam = p.Results.Theta;
            
            if strcmpi(p.Results.MethodName, 'fe')
                A = [0]; b = [1]; s = 1;
            elseif strcmpi(p.Results.MethodName, 'rk4')
                A = [0 0 0 0;
                    1/2 0 0 0;
                    0 1/2 0 0;
                    0 0 1 0]; 
                b = [1/6 1/3 1/3 1/6]; s = 4;
            elseif strcmpi(p.Results.MethodName, 'midpoint')
                A = [0 0;1/2 0]; b = [0 1]; s = 2;
            elseif strcmpi(p.Results.MethodName, 'theta')
                gam = 1/(2*gam);
                A = [0 0;gam 0]; b = [(1 - gam) gam ]; s = 2;
            elseif strcmpi(p.Results.MethodName, 'heuns')
                A = [0 0;1 0]; b = [0.5 0.5]; s = 2;
            elseif strcmpi(p.Results.MethodName, 'ralston')
                A = [0 0; 2/3 0]; b = [1/4 3/4]; s = 2;
            elseif strcmpi(p.Results.MethodName, 'kutta3')
                A = [0 0 0;0.5 0 0;-1 2 0]; b = [1/6 2/3 6/3]; s = 3;
            end
           
            obj = obj@SSPTools.Steppers.ERK(varargin{:},'A',A,'b',b,'s',s);
            obj.name = p.Results.MethodName;
        end
    end
    
end