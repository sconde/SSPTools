classdef IMEXRK < SSPTools.Integrators.RK
    
    properties
       At; bt; ct;
    end
    
    properties ( Access = private)
        isExplicit = false;
        isMSRK = true; 	% Multi-Stage Runge-Kutta
        isButcher = true; % starting with Butcher formulation
        isLowStorage = false;  % need a way to determine is low-storage
        n;
        Y;
        isImplicitLinear = true; %so far handling linear advection
    end
    
    methods
        
        function obj = IMEXRK(varargin)
            obj = obj@SSPTools.Integrators.RK(varargin{:});
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParamValue(p,'name','MSRK-DIRK');
            addParamValue(p, 'isSSP', false);
            addParamValue(p, 'isButcher', true);
            addParamValue(p, 'At', []);
            addParamValue(p, 'bt', []);
            addParamValue(p, 'isLowStorage', false);
            p.parse(varargin{:});

            obj.At = p.Results.At;
            obj.bt = p.Results.bt;
            obj.ct = sum(obj.A,2);
            obj.name = p.Results.name;  
            obj.n = size(obj.y0,1);
            obj.Y = zeros(obj.n, obj.s);
        end
        
        
    end
    
    methods
        function [y] = takeStep(obj, dt)
            
            %check to see if CFL violation
            assert((dt/obj.dx) <= obj.CFL, ...
                sprintf('ERK: CFL Violation (CFL = %3.2f )',dt/obj.dx) );
            
            u0 = obj.y0;
            obj.Y(:,1) = u0;
            
            % intermediate stage value
            for i = 2:obj.s-1
                
                temp = u0;
                for j = 1:i-1
                    temp = temp + dt*obj.A(i,i)*obj.L(dt + obj.c(j), obj.Y(:,j));
                end
                obj.Y(:,i) = temp;
            end
            
            % combine
            y = u0;
            for i = 1:obj.s
                y = y + dt*obj.b(i)*obj.L(dt + obj.c(i), obj.Y(:,i));
            end
            
            obj.y0 = y;
            
        end
        
    end
    
    methods (Access = private)
        
        function  y = linearImplicitStage( y )
        
        end
        
        function y = nonlinearImplicitStage( y )
            
        end
        
    end
    

end
