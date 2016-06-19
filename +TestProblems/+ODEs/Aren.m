classdef Aren < TestProblems.ODEs.ODE
    
    properties
        isLinear;
        mu;
    end
    
    methods
        
        function obj = Aren(varargin)
            
            obj = obj@TestProblems.ODEs.ODE(varargin{:});
            
            p = inputParser;
            p.KeepUnmatched = true;
            
            addParameter(p, 'mu', 0.012277471);
            
            p.parse(varargin{:});
            
            obj.name = 'Aren';
            obj.isLinear = false;
            obj.mu = p.Results.mu;
            
%             %f = zeros(size(y));
%             amu = 0.012277471;
%             amup = 1.0 - amu;
%             sqr = y(1) + amu;
%             r1 = sqr*sqr + y(2)*y(2);
%             r1 = r1 * sqrt(r1);
%             sqr = y(1) - amup;
%             r2 = sqr*sqr + y(2)*y(2);
%             r2 = r2 * sqrt(r2);
            
            obj.f = @(t, u) obj.flux(t, u);
            %obj.f = @obj.flux;
        end
        
    end
    
    methods %(Access = protected)
        
        function f = flux(obj, t, y)
            amup = 1.0 - obj.mu;
            sqr = y(1) + obj.mu;
            r1 = sqr*sqr + y(2)*y(2);
            r1 = r1 * sqrt(r1);
            sqr = y(1) - amup;
            r2 = sqr*sqr + y(2)*y(2);
            r2 = r2 * sqrt(r2);
            
            f = [y(3);
                y(4);
                y(1) + 2.0 * y(4) - amup * (y(1)+obj.mu) / r1 - obj.mu * (y(1)-amup) / r2;
                y(2) - 2.0 * y(3) - amup * y(2) / r1 - obj.mu * y(2) / r2];
        end
        
    end
end
