classdef Mesh < handle
    
    properties
        domain;
        quadratureType
        solutionDegree
        nElements
        nFaces
        nSolutionPoints
    end
    
    properties (Access = protected)
        nNodes;
        solutionPoints
        weights
        elementNodes
        elementFaces
        elementCenter
        elementSize
        nodeCoordinates
        Jacobian
        quadFunction;
        solPts;
    end
    
    methods
        function obj = Mesh(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            % TODO: determine the default values?
            addParameter(p, 'Domain', [-1 1]);
            addParameter(p, 'QuadratureType', []);
            addParameter(p, 'SolutionDegree', []);
            addParameter(p, 'NumberElements', []);
            p.parse(varargin{:});
            
            if ~isempty(p.Results.Domain)
                obj.domain = p.Results.Domain;
            end
            
            if ~isempty(p.Results.QuadratureType)
                obj.quadratureType = p.Results.QuadratureType;
            end
            
            if ~isempty(p.Results.NumberElements)
                obj.nElements = p.Results.NumberElements;
                obj.nFaces    = obj.nElements + 1;
            end
            
            if ~isempty(p.Results.SolutionDegree)
                obj.solutionDegree = p.Results.SolutionDegree;
                obj.nSolutionPoints = obj.solutionDegree + 1;
            end
            
            switch lower(obj.quadratureType)
                case {'legendre'}
                    obj.quadFunction = @(n) obj.GaussLegendre(n);
                case {'lgl'}
                    obj.quadFunction = @(n) obj.GaussLobatto(n);
                otherwise
                    error('quadrature scheme not defined')
            end
        end % constructor
    end
    
    methods  (Access = protected)
        
        % TODO: should probably been moving this to Mesh1D class
        function [x,w] = GaussLegendre(obj, n)
            i  = 1:n-1;
            a  = i./sqrt(4*i.^2-1);
            CM = diag(a,1) + diag(a,-1);
            [V,L] = eig(CM);
            [x,~] = sort(diag(L));
            w = 2 * (V(1,:).^2)';
        end % GaussLegendre
        
        function [x,w,P] = GaussLobatto(obj, kk)
            N = kk-1; % Compute for the number of points kk
            
            % Truncation + 1
            N1=N+1;
            % Use the Chebyshev-Gauss-Lobatto nodes as the first guess
            x=-cos(pi*(0:N)/N)';
            % The Legendre Vandermonde Matrix
            P=zeros(N1,N1);
            % Compute P_(N) using the recursion relation
            % Compute its first and second derivatives and
            % update x using the Newton-Raphson method.
            xold=2;
            while max(abs(x-xold))>eps
                xold=x;
                P(:,1)=1;    P(:,2)=x;
                for k=2:N
                    P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
                end
                x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
            end
            w=2./(N*N1*P(:,N1).^2);
        end % Gausslobatto
       
    end % protected method
    
end
