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
            %TODO: verify all the methods
            
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
                % Hammer & Hollingsworth method
                % Solving Ordinary Differential Equations I:
                % Nonstiff Problems, 2nd Revised Edition
                % E. Hairer, S. P. Norsett, and G. Wanner
                % Table 7.1, pg 205
                A = [0 0;1/3 1/3]; b = [1/4 3/4];
            elseif strcmpi(p.Results.MethodName, 'sdirk22')
                % Computer Methods for ODEs and DAEs
                % U. M. Ascher and L. R. Petzold
                % p. 106\n
                % gamma = (2+-sqrt(2))/2\n
                gam = (2 + sqrt(2))/2;
                A = [gam 0;1-gam gam]; b = [1-gam gam];
            elseif strcmpi(p.Results.MethodName, 'sdirk23')
                % Solving Ordinary Differential Equations I:
                % Nonstiff Problems
                % 2nd Revised Edition
                % E. Hairer and G. Wanner
                % Table 7.2, pg 207
                gam = (3 + sqrt(3))/6;
                A = [gam 0; (1-2*gam) gam]; b = [1/2 1/2];
            elseif strcmpi(p.Results.MethodName, 'sdirk55')
                % Solving Ordinary Differential Equations II:
                % Stiff and Differential-Algebraic Problems
                % 2nd Revised Edition
                % E. Hairer and G. Wanner
                % pg100
                A = zeros(5); gam = (6 - sqrt(6))/10;
                A(2,1) = (-6 + 5*sqrt(6))/14;
                A(3,1) = (888 + 607*sqrt(6))/2850;
                A(3,2) = (126 - 161*sqrt(6))/1425;
                A(4,1) = (3153 - 2082*sqrt(6))/14250;
                A(4,2) = (3213 + 1148*sqrt(6))/28500;
                A(4,3) = (-267 + 88*sqrt(6))/500;
                A(5,1) = (-32583 + 14838*sqrt(6))/71250;
                A(5,2) = (-17199 + 364*sqrt(6))/142500;
                A(5,3) = (1329 - 544*sqrt(6))/2500;
                A(5,4) = (-96 + 131*sqrt(6))/625;
                A = A + gam*eye(5);
                b = [0 0 1/9 (16 - sqrt(6))/36 (16+sqrt(6))/36];
            elseif strcmpi(p.Results.MethodName, 'sdirk54')
                % Solving Ordinary Differential Equations II:
                % Stiff and Differential-Algebraic Problems
                % 2nd Revised Edition
                % E. Hairer and G. Wanner
                % pg100
                A = zeros(5); gam = 1/4;
                A(2,1) = 1/2;
                A(3,1) = 17/50;
                A(3,2) = -1/25;
                A(4,1) = 371/1360;
                A(4,2) = -137/2720;
                A(4,3) = 15/544;
                A(5,1) = 25/24;
                A(5,2) = -49/48;
                A(5,3) = 125/16;
                A(5,4) = -85/12;
                A = A + gam*eye(5);
                b = [25/24 -49/28 125/16 -85/12 1/4];
            elseif strcmpi(p.Results.MethodName, 'sdirk34')
                % A-stable
                % Solving Ordinary Differential Equations II:
                % Stiff and Differential-Algebraic Problems
                % 2nd Revised Edition
                % E. Hairer and G. Wanner
                % pg100
                % gamma = (1/sqrt(3))*cos(pi/18)+1/2\n"
                % delta = 1/(6*(2*gamma-1)^2)\n"
                gam = (1/sqrt(3))*cos(pi/18)+1/2;
                delt = 1/(6*(2*gam - 1)^2);
                A = zeros(3);
                A(2,1) = 1/2 - gam;
                A(3,1) = 2*gam;
                A(3,2) = 1 - 4*gam;
                A = A + gam*eye(3);
                b = [delt (1-2*delt) delt];
            end
            
            obj = obj@SSPTools.Steppers.DIRK(varargin{:},'A',A,'b',b);
            obj.name = p.Results.MethodName;
        end
    end
    
end
