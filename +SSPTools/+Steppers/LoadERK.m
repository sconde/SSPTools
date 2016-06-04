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
                A = [0]; b = [1]; r = 1;
            elseif strcmpi(p.Results.MethodName, 'ssp22')
                A = [0 0;1 0]; b = [1/2 1/2]; r = 1;
            elseif strcmpi(p.Results.MethodName, 'ssp33')
                A = [0 0 0;1 0 0;1/4 1/4 0]; b = [1/6 1/6 2/3]; r = 1;
            elseif strcmpi(p.Results.MethodName, 'ssp54')
                A = zeros(5); r = 1.508180049189838;
                A(tril(true(5),-1)) = [0.391752226571889; 0.217669096261169; ...
                    0.082692086657811; 0.067966283637115; 0.368410593050372; ...
                    0.139958502191896; 0.115034698504632; 0.251891774271693; ...
                    0.207034898597385; 0.544974750228520];
                b = [0.146811876084786; 0.248482909444976; ...
                     0.104258830331980; 0.274438900901350; 0.226007483236908]; 
            elseif strcmpi(p.Results.MethodName, 'rk4')
                A = [0 0 0 0;
                    1/2 0 0 0;
                    0 1/2 0 0;
                    0 0 1 0]; 
                b = [1/6 1/3 1/3 1/6]; r = -Inf;
            elseif strcmpi(p.Results.MethodName, 'midpoint')
                A = [0 0;1/2 0]; b = [0 1]; r = -Inf;
            elseif strcmpi(p.Results.MethodName, 'theta')
                gam = 1/(2*gam);
                A = [0 0;gam 0]; b = [(1 - gam) gam ]; r = -Inf;
            elseif strcmpi(p.Results.MethodName, 'heuns')
                A = [0 0;1 0]; b = [0.5 0.5]; r = -Inf;
            elseif strcmpi(p.Results.MethodName, 'ralston')
                A = [0 0; 2/3 0]; b = [1/4 3/4]; r = -Inf;
            elseif strcmpi(p.Results.MethodName, 'kutta3')
                A = [0 0 0;0.5 0 0;-1 2 0]; b = [1/6 2/3 6/3]; 
                r = -Inf;
            end
           
            obj = obj@SSPTools.Steppers.ERK(varargin{:},'A',A,'b',b,'r',r);
            obj.name = p.Results.MethodName;
        end
    end
    
end