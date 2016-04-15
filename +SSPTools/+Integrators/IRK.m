classdef IRK < SSPTools.Integrators.RK
    
    properties
        A = []; b = []; c = []; alpha = []; s = []; r = [];
        isExplicit = true;
        isMSRK = true; 	% Multi-Stage Runge-Kutta
        isButcher = true; % starting with Butcher formulation
        isSSP = false;
        isLowStorage = false;  % need a way to determine is low-storage
        steps = 1; % single-step methods
        n;
    end
    
    methods
        
        function obj = IRK(varargin)
            obj = obj@SSPTools.Integrators.RK(varargin{:});
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParamValue(p,'name','MSRK-ERK');
            addParamValue(p, 'A', []);
            addParamValue(p, 'b', []);
            addParamValue(p, 's', []);
            addParamValue(p, 'r', []);
            addParamValue(p, 'alpha', []);
            addParamValue(p, 'isSSP', false);
            addParamValue(p, 'isButcher', true);
            addParamValue(p, 'isLowStorage', false);
            addParamValue(p, 'dydt', []);
            addParamValue(p, 'y0', []);
            p.parse(varargin{:});
                        
            assert(isequal(p.Results.s, size(p.Results.A,1)),...
                sprintf('ERK A:Stage-count -- Num-Rows(A) != %d',p.Results.s));
            
            obj.A = p.Results.A;
            obj.b = p.Results.b;
            obj.c = sum(obj.A,2);
            obj.s = p.Results.s;
            obj.alpha = p.Results.alpha;
            obj.name = p.Results.name;
            
            if isa(p.Results.dydt, 'function_handle')
                obj.dydt = p.Results.dydt;
            end
            
            if ~isempty(p.Results.y0)
                obj.y0 = p.Results.y0(:);
                obj.n = size(obj.y0,1);
            end
                
        end
        
        
    end
    
    methods %( Access = protected )
        function [y] = takeStep(obj, dt)
            
            u0 = obj.y0;
            Y = zeros(obj.s, obj.n);
            Y(1) = u0;
            
            % intermediate stage value
            for i = 2:obj.s-1
                
                temp = u0;
                for j = 1:i-1
                    temp = temp + dt*obj.A(i,i)*obj.dydt(dt + obj.c(j), Y(j));
                end
                Y(i) = temp;
            end
            
            % combine
            y = u0;
            for i = 1:obj.s
                y = y + dt*obj.b(i)*obj.dydt(dt+obj.c(i), Y(i));
            end
            
            obj.y0 = y;
            
        end
        
    end
    
    methods ( Access = private )
        
        
    end
end
