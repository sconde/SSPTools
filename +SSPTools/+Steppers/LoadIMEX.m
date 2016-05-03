classdef LoadIMEX < SSPTools.Steppers.IMEXRK
    
    properties
        MethodName;
    end
    
    properties (Access = private)
        gam;
    end
    
    methods
        function obj = LoadIMEX(varargin)
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'MethodName','SSP11LPM');
            addParameter(p, 'Gamma', 0.1);
            p.parse(varargin{:});
            
            gam = p.Results.Gamma;
            
            if strcmpi(p.Results.MethodName, 'ssp11lpm')
                A = [0]; b = [1]; At = [1]; bt = [1]; s = 1;
            elseif strcmpi(p.Results.MethodName, 'arslpum')
                A = [0 0;1 0]; b = [1 0]; s = 2;
                At = [0 0;0 1]; bt = [0 1];
            elseif strcmpi(p.Results.MethodName, 'stormer-verlet')
                A = [0 0;0.5 0.5]; b = [0.5 0.5]; s = 2;
                At = [0.5 0;0.5 0]; bt = b;
            elseif strcmpi(p.Results.MethodName, 'ssp22')
                A = [0 0;1 0]; b = [0.5 0.5]; s = 2;
                At = [gam 0; 1-2*gam gam]; bt = b;
            elseif strcmpi(p.Results.MethodName, 'ssp22um')
                A = [0 0;1 0]; b = [0.5 0.5]; s = 2;
                At = [0 0; 0.5 0.5]; bt = b;
            elseif strcmpi(p.Results.MethodName, 'ssp22lum')
                A = [0 0 0; 0.5 0 0;0.5 0.5 0]; b = [1/3 1/3 1/3]; s = 3;
                At = [1/5 0 0; 1/10 1/5 0;1/3 1/3 1/3]; bt = b;
            elseif strcmpi(p.Results.MethodName, 'ssp22stiff')
                A = [0 0 0;0.5 0 0;0.5 0.5 0]; b = [1/3 1/3 1/3]; s = 3;
                At = [1/4 0 0;0 1/4 0;1/3 1/3 1/3]; bt = b;
            elseif strcmpi(p.Results.MethodName, 'ssp22lstable')
                A = [0 0 0;1 0 0;1/4 1/4 0]; b = [1/6 1/6 2/3]; s = 3;
                At = [gam 0 0;1-2*gam gam 0 0.5-gam 0 gam];
                bt = b;
            elseif strcmpi(p.Results.MethodName, 'peaceman-rachford')
                A = [0 0 0;0.5 0 0;0.5 0 0.5]; b = [0.5 0 0.5]; s = 3;
                At = [0 0 0; 0 0.5 0;0 1 0]; bt = [0 1 0];
            elseif strcmpi(p.Results.MethodName, 'ars')
                A = [0 0 0 0 0;
                    1/2 0 0 0 0;
                    11/18 1/18 0 0 0;
                    5/6 -5/6 1/2 0 0;
                    1/4 7/4 3/4 -7/4 0];
                b = [1/4 7/4 3/4 -7/4 0]; s = 5;
                At = [0 0 0 0 0;
                    0 1/2 0 0 0;
                    0 1/6 1/2 0 0;
                    0 -1/2 1/2 1/2 0;
                    0 3/2 -3/2 1/2 1/2];
                bt = [0 3/2 -3/2 1/2 1/2];
            elseif strcmpi(p.Results.MethodName', 'bpr')
                A = [0 0 0 0 0;
                    1 0 0 0 0;
                    4/9 2/9 0 0 0;
                    1/4 0 3/4 0 0;
                    1/4 0 3/4 0 0];
                b = [1/4 0 3/4 0 0]; s = 5;
                At = [0 0 0 0 0;
                    1/2 1/2 0 0 0;
                    5/128 -1/9 1/2 0 0;
                    1/2 0 0 1/2 0;
                    1/4 0 3/4 -1/2 1/2];
                bt = [1/4 0 3/4 -1/2 -1/2];
            elseif strcmpi(p.Results.MethodName, 'imex1')
                A = [0 0;0.5 0]; b = [0 1]; s = 2;
                At = [0 0;0 0.5]; bt = [0 0.5];
            end

            obj = obj@SSPTools.Steppers.IMEXRK(varargin{:},'A',A,'b',b,'s',s,'At',At,'bt',bt);
            obj.name = p.Results.MethodName;
        end
    end
    
end