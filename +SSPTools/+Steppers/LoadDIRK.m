classdef LoadDIRK < SSPTools.Steppers.DIRK
    %FIXME: the methods are explicit
    
    properties
        MethodName;
    end
    
    properties (Access = private)
        gam;
    end
    
    methods
        function obj = LoadDIRK(varargin)
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'MethodName','FE');
            addParameter(p, 'Constant', 0.1);
            p.parse(varargin{:});
            
            if ~isempty(p.Results.Constant)
                gam = p.Results.Constant;
            end
            
            if strcmpi(p.Results.MethodName, 'be')
                A = [1]; b = [1]; s = 1;
            elseif strcmpi(p.Results.MethodName, 'trapezoid')
                A = [0 0;1/2 1/2]; b = [1/2 1/2]; s = 2;
            elseif strcmpi(p.Results.MethodName, 'dirk23')
                A = [0 0;1/3 1/3]; b = [1/4 3/4];
            elseif strcmpi(p.Results.MethodName, 'sdirk22')
                gam = (2 + sqrt(2))/2;
                A = [gam 0;1-gam gam]; b = [1-gam gam];
            elseif strcmpi(p.Results.MethodName, 'sdirk23')
                gam = (3 + sqrt(3))/6;
                A = [gam 0; (1-2*gam) gam]; b = [1/2 1/2];
            end
            
            obj = obj@SSPTools.Steppers.DIRK(varargin{:},'A',A,'b',b);
            obj.name = p.Results.MethodName;
        end
    end
    
end