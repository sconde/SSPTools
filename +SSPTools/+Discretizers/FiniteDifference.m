classdef FiniteDifference < SSPTools.Discretizers.Discretize
    
    properties
        isSSP = true;
        name = 'Finite Difference';
        bc;
        nx;
        domain;
        D;
        dx;
        systemSize;
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
            addParameter(p, 'domain', []);
            addParameter(p, 'Problem', []);
            p.parse(varargin{:});
            
            if ~isempty(p.Results.Problem)
                obj.problem = p.Results.Problem;
            end
            
            obj.isSSP = p.Results.isSSP;
            obj.bc = p.Results.bc;
            obj.derivativeOrder = p.Results.derivativeOrder;
            
            if ~obj.problem.isSystem
                
                obj.domain = p.Results.domain;
                obj.nx = p.Results.N;
                
                if ~isempty(obj.domain)
                    obj.x = linspace(obj.domain(1),obj.domain(2),obj.nx);
                    obj.dx = obj.x(2) - obj.x(1);
                end
                
                obj.systemSize = 1;
            else
                obj.dx      = obj.problem.dx;
                obj.x       = obj.problem.x;
                obj.nx       = obj.problem.N;
                obj.domain  = obj.problem.domain;
                obj.systemSize = obj.problem.systemSize;
            end
            
            obj.nx = numel(obj.x);
            obj.x = obj.x(:);
            % define the descritization matrix
            [r,c] = deal(zeros(obj.nx,1));
            r([1 obj.nx]) = [1 -1];
            c([1 2]) = [1 -1];
            obj.D = toeplitz(c,r)/(-obj.dx);
                        
        end % end constructor
        
    end
    
    methods
        function [y] = L(obj, t, y)
            yx = zeros(size(y));
            IDX = obj.nx*[0:obj.systemSize-1 ; 1:obj.systemSize]';
            IDX = IDX + [ones(size(IDX(:,1))) zeros(size(IDX(:,1)))];
            fx = obj.problem.f(t,y);
            for i = 1:obj.systemSize
                idx = IDX(i,1); idx_ub = IDX(i,2);
                yx(idx:idx_ub) = obj.D*fx(idx:idx_ub);
            end
            y = yx;            
        end
    end
    
    
end
