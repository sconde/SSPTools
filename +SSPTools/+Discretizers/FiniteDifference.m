classdef FiniteDifference < SSPTools.Discretizers.Discretize
    
    properties
        isSSP = true;
        name = 'Finite Difference';
        bc;
        nx;
        domain;
        D;
        dx;
    end
    
    methods
        
        function obj = FiniteDifference(varargin)
            obj = obj@SSPTools.Discretizers.Discretize();
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'isSSP', true);
            addParameter(p, 'derivativeOrder', 1);
            addParameter(p, 'bc', 'periodic');
            addParameter(p, 'N', 10);
            addParameter(p, 'domain', [-1 1]);
            addParameter(p, 'Problem', []);
            p.parse(varargin{:});
            
            obj.isSSP = p.Results.isSSP;
            obj.bc = p.Results.bc;
            obj.derivativeOrder = p.Results.derivativeOrder;
            obj.nx = p.Results.N;
            obj.domain = p.Results.domain;
            
            if ~isempty(p.Results.Problem)
                obj.problem = p.Results.Problem;
            end
            
            if ~isempty(obj.domain)
                obj.x = linspace(obj.domain(1),obj.domain(2),obj.nx);
                obj.dx = obj.x(2) - obj.x(1);
                obj.nx = numel(obj.x);
                obj.x = obj.x(:);
            end
            
            % define the descritization matrix
            [r,c] = deal(zeros(obj.nx,1));
            r([1 obj.nx]) = [1 -1];
            c([1 2]) = [1 -1];
            obj.D = toeplitz(c,r)/(-obj.dx);
            
        end
        
    end
    
    methods 
        function [y] = L(obj, y)
            y = obj.D*y;
        end
    end
    
    
end
