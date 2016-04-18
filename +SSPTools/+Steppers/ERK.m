classdef ERK < SSPTools.Steppers.RK
    
    properties
        
    end
    
    properties ( Access = private)
        isExplicit = true;
        isMSRK = true; 	% Multi-Stage Runge-Kutta
        isButcher = true; % starting with Butcher formulation
        isLowStorage = false;  % need a way to determine is low-storage
        n;
        Y;
        u0;
    end
    
    methods
        
        function obj = ERK(varargin)
            obj = obj@SSPTools.Steppers.RK(varargin{:});
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParamValue(p,'name','MSRK-ERK');
            addParamValue(p, 'isSSP', false);
            addParamValue(p, 'isButcher', true);
            addParamValue(p, 'isLowStorage', false);
            p.parse(varargin{:});

            obj.name = p.Results.name;  
            obj.n = size(obj.x,1);
            obj.Y = zeros(obj.n, obj.s);
            if isa(obj.y0, 'function_handle')
                obj.u0 = obj.y0(obj.x);
            else
                obj.u0 = obj.y0;
            end
            
        end
        
        
    end
    
    methods %( Access = protected )
        function [y] = takeStep(obj, dt)
            
            %check to see if CFL violation
            assert((dt/obj.dx) <= obj.CFL, ...
                sprintf('ERK: CFL Violation (CFL = %3.2f )',dt/obj.dx) );
            
            u0 = obj.u0;
            obj.Y(:,1) = u0;
            
            % intermediate stage value
            for i = 2:obj.s-1
                
                temp = u0;
                for j = 1:i-1
                    temp = temp + dt*obj.A(i,j)*obj.L(dt + obj.c(j), obj.Y(:,j));
                end
                obj.Y(:,i) = temp;
            end
            
            % combine
            y = u0;
            for i = 1:obj.s
                y = y + dt*obj.b(i)*obj.L(dt + obj.c(i), obj.Y(:,i));
            end
            
            obj.u0 = y;
            
        end
        
    end
    

end
