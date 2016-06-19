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
                %rk = load('dormandprince54.mat');
                A = zeros(7);
                A(2,1)=0.2;
                A(3,1)=3.0/40.0;
                A(3,2)=9.0/40.0;
                A(4,1)=44.0/45.0;
                A(4,2)=-56.0/15.0;
                A(4,3)=32.0/9.0;
                A(5,1)=19372.0/6561.0;
                A(5,2)=-25360.0/2187.0;
                A(5,3)=64448.0/6561.0;
                A(5,4)=-212.0/729.0;
                A(6,1)=9017.0/3168.0;
                A(6,2)=-355.0/33.0;
                A(6,3)=46732.0/5247.0;
                A(6,4)=49.0/176.0;
                A(6,5)=-5103.0/18656.0;
                A(7,1)=35.0/384.0;
                A(7,3)=500.0/1113.0;
                A(7,4)=125.0/192.0;
                A(7,5)=-2187.0/6784.0;
                A(7,6)=11.0/84.0;
                
                b = zeros(7,1);
                b(1) = 35/384;
                b(3) = 500/1113;
                b(4) = 125/192;
                b(5) = -2187/6784;
                b(6) = 11/84;
                
                bhat = zeros(7,1);
                bhat(1) = 5179/57600;
                bhat(3) = 7571/16695;
                bhat(4) = 393/640;
                bhat(5) = -92097/339200;
                bhat(6) = 187/2100;
                bhat(7) = 1/40;
                
                c = sum(A,2);
                
                E = zeros(size(c));
                E(1)=71.0/57600.0;
                E(3)=-71.0/16695.0;
                E(4)=71.0/1920.0;
                E(5)=-17253.0/339200.0;
                E(6)=22.0/525.0;
                E(7)=-1.0/40.0;
                
                D = zeros(size(E));
                D(1)=-12715105075.0/11282082432.0;
                D(3)=87487479700.0/32700410799.0;
                D(4)=-10690763975.0/1880347072.0;
                D(5)=701980252875.0/199316789632.0;
                D(6)=-1453857185.0/822651844.0;
                D(7)=69997945.0/29380423.0;
                
                %A = rk.A; b = rk.b(:); bhat = rk.bhat(:); 
                p = 4; phat = 5;
                
            elseif strcmpi(inpPar.Results.MethodName, '38rule43')
                % 3/8-rule fourth-order method
                % equipped it with the embedded formula (4.9) of order 3
                % Solving Ordinary Differential Equation I
                % Hairer. pg. 170 (Fig. 4.1)
                p = 4; phat = 3;
                A = zeros(5);
                A(2,1) = 1/3;
                A(3,1) = -1/3; A(3,2) = 1;
                A(4,1) = 1; A(4,2) = -1; A(4,3) = 1;
                A(5,1) = 1/8; A(5,2) = 3/8; A(5,3) = 3/8; A(5,4) = 1/8;
                b = A(5,:);
                bhat = [1/12 1/2 1/4 0 1/6];
            end
            
            obj = obj@SSPTools.Steppers.EmbeddedERK(varargin{:},'A',A,'b',b,'r',r,...
                'bhat',bhat, 'p', p, 'phat', phat);
            obj.name = inpPar.Results.MethodName;
        end
    end
    
end
