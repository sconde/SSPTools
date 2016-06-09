classdef LoadEmbeddedERK < SSPTools.Steppers.EmbeddedERK
    
    properties
        MethodName;
    end
    
    properties (Access = private)
        gam;
    end
    
    methods
        function obj = LoadEmbeddedERK(varargin)
            
            inpPar = inputParser;
            inpPar.KeepUnmatched = true;
            addParameter(inpPar,'MethodName','Merson45');
            inpPar.parse(varargin{:});
            
            
            r = [];
            %TODO: where is the 38Ruler???
            
            if strcmpi(inpPar.Results.MethodName, 'merson45')
                % Solving Ordinary Differential Equation I
                % Hairer, pg167. table 4.1
                p = 4; phat = 5;
                A = [0 0 0 0 0;
                    1/3 0 0 0 0;
                    1/6 1/6 0 0 0;
                    1/8 0 3/8 0 0;
                    1/2 0 -3/2 2 0];
                b = [1/6 0  0 2/3 1/6];
                bhat = [1/10 0 3/10 2/5 1/5];
            elseif strcmpi(inpPar.Results.MethodName, 'zonneveld43')
                % Solving Ordinary Differential Equation I
                % Hairer, pg167. table 4.1
                p = 4; phat = 3;
                A = [0 0 0 0 0;
                    1/2 0 0 0 0;
                    0 1/2 0 0 0;
                    0 0 1 0 0;
                    5/32 7/32 12/32 -1/32 0];
                b = [1/6 1/3 1/3 1/6 0];
                bhat = [-1/2 7/3 7/3 13/6 -16/3];
            elseif strcmpi(inpPar.Results.MethodName, 'felhberg45')
                % https://en.wikipedia.org/wiki/Runge-Kutta-Fehlberg_method
                p = 4; phat = 5;
                A = zeros(6);
                A(2,1) = 1/4;
                A(3,1) = 3/32; A(3,2) = 9/32;
                A(4,1) = 1932/2197; A(4,2) = -7200/2197; A(4,3) = 7296/2197;
                A(5,1) = 439/216; A(5,2) = -8; A(5,3) = 3680/513; A(5,4) = -845/4104;
                A(6,1) = -8/27; A(6,2) = 2; A(6,3) = -3544/2565; A(6,4) = 1859/4104;
                A(6,5) = -11/40;
                bhat = [16/135 0 6656/12825 28561/56430 -9/50 2/55];
                b = [25/216 0 1408/2565 2197/4104 -1/5 0];
            elseif strcmpi(inpPar.Results.MethodName, 'dormandprince54')
                %TODO: confirm order of accuracy
                % https://en.wikipedia.org/wiki/Dormand-Prince_method
                A = zeros(7);
                A(2,1) = 1/5;
                A(3,1) = 3/40; A(3,2) = 9/40;
                A(4,1) = 44/35; A(4,2) = -56/15; A(4,3) = 32/9;
                A(5,1) = 19372/6561; A(5,2) = -25360/2187; A(5,3) = 64448/6561;
                A(5,4) = -212/729;
                A(6,1) = 9017/3168; A(6,2) = -355/33; A(6,3) = 46732/5247;
                A(6,4) = 49/176; A(6,5) = -5103/18656;
                A(7,1) = 35/384; A(7,2) = 0; A(7,3) = 500/1113; A(7,4) = 125/192;
                A(7,5) = -2187/6784; A(7,6) = 11/84;
                b = A(7,:);
                bhat = [5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40];
            end
           
            obj = obj@SSPTools.Steppers.EmbeddedERK(varargin{:},'A',A,'b',b,'r',r,...
                'bhat',bhat, 'p', p, 'phat', phat);
            obj.name = inpPar.Results.MethodName;
        end
    end
    
end
