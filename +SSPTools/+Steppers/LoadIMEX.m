classdef LoadIMEX < SSPTools.Steppers.IMEXRK
    % kupka1 :Total-Variation-Diminishing Implicit Explicit Runge-Kutta
    % Methods for the Simulation of Double-Diffusive Convection in
    % Astrophysics. Kupka, Happenhofer, Higueras, Koch. 2011.
    
    % higuera1: 
    
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
            addParameter(p,'MethodName','SSP1111LPM');
            addParameter(p, 'Gamma', 0.1);
            p.parse(varargin{:});
            
            gam = p.Results.Gamma;
            r = 0; rt = 0;
                        
            if strcmpi(p.Results.MethodName, 'imexssp1111lpm')
                %\cite[higuera1, eq. 9.1]
                A = [0]; b = [1]; At = [1]; bt = [1]; s = 1;
                r = 1; rt = Inf;
            elseif strcmpi(p.Results.MethodName, 'arslpum')
                A = [0 0;1 0]; b = [1 0]; s = 2;
                At = [0 0;0 1]; bt = [0 1];
            elseif strcmpi(p.Results.MethodName, 'stormer-verlet')
                A = [0 0;0.5 0.5]; b = [0.5 0.5]; s = 2;
                At = [0.5 0;0.5 0]; bt = b;
            elseif strcmpi(p.Results.MethodName, 'imexssp2222')
                % \cite[kupka1, higueras1]
                gam = 1 - 1/sqrt(2);
                A = [0 0;1 0]; b = [0.5 0.5]; s = 2;
                At = [gam 0; 1-2*gam gam]; bt = b;
            elseif strcmpi(p.Results.MethodName, 'imexssp2222pm')
                % \cite[kupka1, higueras1]
                gam = 0.24;
                A = [0 0;1 0]; b = [0.5 0.5]; s = 2;
                At = [gam 0; 1-2*gam gam]; bt = b;
                r = 1; rt = 3.57;                
            elseif strcmpi(p.Results.MethodName, 'ssp22um')
                A = [0 0;1 0]; b = [0.5 0.5]; s = 2;
                At = [0 0; 0.5 0.5]; bt = b;
            elseif strcmpi(p.Results.MethodName, 'imexssp2332lum')
                % \cite[kupka1,higueras1 eq. 9.3]
                A = [0 0 0; 0.5 0 0;0.5 0.5 0]; b = [1/3 1/3 1/3]; s = 3;
                At = [1/5 0 0; 1/10 1/5 0;1/3 1/3 1/3]; bt = b;
                r = 2; rt = (5/9)*(sqrt(70) - 4);
            elseif strcmpi(p.Results.MethodName, 'imexssp2332lspum')
                % \cite[higueras eq. 3.30]
                A = [0 0 0;5/6 0 0;11/24 11/24 0]; s = 3;
                b = [24/25 1/5 4/11]; bt = b;
                At = [2/11 0 0;205/462 2/11 0;2033/4620 21/110 2/11];
                r = 1.2; rt = 3.82;
            elseif strcmpi(p.Results.MethodName, 'imexssp2332lpum')
                % \cite[higueras eq. 4.2]
                A = [0 0 0;1/2 0 0;1/2 1/2 0]; s = 3;
                b = [1/3 1/3 1/3]; bt = b;
                At = [2/11 0 0;41/154 2/11 0;289/847 42/121 2/11];
                r = 2; rt = 3.09;
            elseif strcmpi(p.Results.MethodName, 'imexssp2332lpm')
                % \cite[higueras eq. 4.5]
                A = [0 0 0;1/2 0 0;1/2 1/2 0]; s = 3;
                b = [1/3 1/3 1/3]; bt = b;
                At = [2/11 0 0;2583/13310 2/11 0;39731/139755 10/21 2/11];
                r = 2; rt = 2.34;
            elseif strcmpi(p.Results.MethodName, 'imexssp2332lu')
                % \cite[higueras eq. 9.4]
                A = [0 0 0;1/2 0 0;1/2 1/2 0]; s = 3;
                b = [1/3 1/3 1/3]; bt = b;
                At = [1/4 0 0;0 1/4 0;1/3 1/3 1/3];
                r = 2; rt = 2.4;
            elseif strcmpi(p.Results.MethodName, 'ssp22stiff')
                A = [0 0 0;0.5 0 0;0.5 0.5 0]; b = [1/3 1/3 1/3]; s = 3;
                At = [1/4 0 0;0 1/4 0;1/3 1/3 1/3]; bt = b;
            elseif strcmpi(p.Results.MethodName, 'ssp22lstable')
                A = [0 0 0;1 0 0;1/4 1/4 0]; b = [1/6 1/6 2/3]; s = 3;
                At = [gam 0 0;1-2*gam gam 0 0.5-gam 0 gam];
                bt = b;
            elseif strcmpi(p.Results.MethodName, 'imexssp2332lu')
                % \cite[higueras1]
                A = [0 0 0;1/2 0 0;1/2 1/2 0;1/3 1/3 1/3]; s = 3;
                b = [1/3 1/3 1/3];
                At = [1/4 0 0;0 1/4 0;1/3 1/3 1/3]; bt = b;
            elseif strcmpi(p.Results.MethodName, 'imexssp3333')
                % \cite[kupka1]
                A = [0 0 0;1 0 0;1/4 1/4 0]; b = [1/6 1/6 2/3]; s = 3;
                At = [0 0 0;14/15 1/15 0;7/30 1/5 1/15]; bt = [1/6 1/6 2/3];
                r = 1; rt = (5/47)*(13 - 2*sqrt(7));
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
                % TODO: fix this, somehow not catching this
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

            obj = obj@SSPTools.Steppers.IMEXRK(varargin{:},'A',A,'b',b,...
                's',s,'At',At,'bt',bt,'r',r, 'rt',rt);
            obj.name = p.Results.MethodName;
        end
    end
    
end