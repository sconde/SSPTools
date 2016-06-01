classdef NDGBase < handle
    
    properties
        kDeg
        nNodes
        x;
        mesh
        problem
        dx;
    end
    
    properties (Access = protected)
        Vandermonde
        Vandermonde2
        GradVandermonde
        GradVandermonde2
        MassMatrix
        invMassMatrix
        CoefDiffMatrix
        legLeftEnd
        legRightEnd
        nodalMassMatrix
        nodalInvMassMatrix
        nodalCoefDiffMatrix
    end
    
    methods
        function obj = NDGBase(varargin)
            % TODO: determine the default values?
            p = inputParser;
            p.KeepUnmatched = true;
            
            addParameter(p, 'Mesh', []);
            addParameter(p, 'Problem', []);
            p.parse(varargin{:});
            
            assert(~isempty(p.Results.Mesh),'Need a Mesh to work with');
            
            if isa(p.Results.Mesh, 'NDG.Mesh1D')
                obj.mesh = p.Results.Mesh;
            else
                error('Only 1D Mesh is supported right now');
            end
            
            if ~isempty(p.Results.Problem)
                obj.problem = p.Results.Problem;
            end
            
            obj.nNodes = length(obj.mesh.solPts);
            obj.kDeg = length(obj.mesh.solPts)-1;
            obj.x = obj.mesh.solPts;
            
            % construct the Vandermond Matrix
            V = zeros(obj.kDeg);% Allocate
            for l = 0:obj.kDeg 	% All Polynomial Degrees up to kDeg
                j = l + 1;      % Dummy index
                for i = 1:obj.nNodes  % All Element Points in xi
                    V(i,j) = obj.legendreP(obj.x(i),l);
                end
            end
            obj.Vandermonde = V;
            
            % Construct Legendre Vandermonde Matrix, V,
            V = zeros(obj.kDeg);% Allocate
            for l = 0:obj.kDeg 	% All Polynomial Degrees up to kDeg
                j = l + 1;      % Dummy index
                for i = 1:obj.nNodes  % All Element Points in xi
                    V(i,j) = obj.JacobiP(obj.x(i),0,0,l);
                end
            end
            obj.Vandermonde2 = V;
            
            l = (0:obj.kDeg);
            obj.legLeftEnd = (-1).^l;
            obj.legRightEnd = 1.^l;
            obj.MassMatrix = diag(2./(2*l+1));
            obj.invMassMatrix = inv(obj.MassMatrix);
            
            
            k = obj.kDeg; Dhat = zeros(k+1);
            for m = 0:k
                for j = m+1:2:k
                    Dhat(m+1,j+1) = (2*m+1);
                end
            end
            
            obj.CoefDiffMatrix = Dhat;
            
            %function nM = get.nodalMassMatrix(obj)
            obj.nodalMassMatrix = inv(obj.Vandermonde2*obj.Vandermonde2');
            %end
            
            %function ninvM = get.nodalInvMassMatrix(obj)
            obj.nodalInvMassMatrix = inv(obj.nodalMassMatrix);
            %end
            
            %function Vr = get.GradVandermonde(obj)
            Vr = zeros(obj.kDeg);% Allocate
            for l = 0:obj.kDeg 	% All Polynomial Degrees up to kDeg
                j = l + 1;      % Dummy index
                for i = 1:obj.nNodes  % All Element Points in xi
                    Vr(i,j) = obj.dJacobiP(obj.x(i),0,0,l);
                end
            end
            obj.GradVandermonde = Vr;
            %end
            
            %function nDr = get.nodalCoefDiffMatrix(obj)
            obj.nodalCoefDiffMatrix = obj.GradVandermonde/obj.Vandermonde2;
            %end
            
        end % NDG constructor
        
    end
    
    methods (Access = protected)
        
        function dF = residual(obj, t, u, flux,dflux,Lift,Dr)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %       Compute the Residual for 1d wave equation using Nodal DG
            %
            %                       Residual = dF/dxi
            %                 where F = is our Flux function
            %
            %              coded by Manuel Diaz, NTU, 2013.10.29
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % compute fluxes in node coordinates
            f = flux(t, u);
            
            % Interpolate u and flux values at the boundaries of Ij
            % use only: {'LGL','ChebyshevMod'}
            u_lbd = u(1,:);
            u_rbd = u(end,:);
            f_lbd = f(1,:);
            f_rbd = f(end,:);
            
            % Build Numerical fluxes across faces
            u_pface = [u_lbd,0]; % + side
            u_nface = [0,u_rbd]; % - side
            
            % Apply Periodic BCs
            u_nface(1) = u_nface(end); % left BD
            u_pface(end) = u_pface(1); % right BD
            
            % Apply Neumann BCs
            %u_nface(1) = u_pface(1); % left BD
            %u_pface(end) = u_nface(end);% u_pface(end); % right BD
            
            % LF numerical flux
            alpha = max(max(abs(dflux(u))));
            nflux = 0.5*(flux(t, u_nface)+flux(t, u_pface)-alpha*(u_pface-u_nface));
            nfluxL = nflux(1:end-1); nfluxR = nflux(2:end);
            
            % Compute the derivate: F = f + lL*(nfluxL-f_bdL) - lR*(nfluxR-f_bdR)
            dF = -Dr*f + Lift*[(nfluxL-f_lbd);(f_rbd-nfluxR)];
            
            dF = -dF;
        end
        
        function legP = legendreP(obj, xi,l)
            % Construct Legendre Polynomials
            %**************************************************************
            % Compute the value of the legendre polinomial of degree 'kDeg'
            % for every column value 'xi'
            %**************************************************************
            x = sym('x');
            temp = simplify(1/(2^l*factorial(l))*diff((x^2-1)^l,x,l));
            legP = subs(temp,xi);
        end
        
        function dlegP = dlegendreP(obj, xi,l)
            % Construct derivatives of Legendre Polynomials
            %**************************************************************
            % Compute the derivative value of the legendre polinomial of
            % degree 'kDeg' for every column value 'xi'
            %**************************************************************
            x = sym('x');
            legP = simplify(1/(2^l*factorial(l))*diff((x^2-1)^l,x,l));
            dlegP = subs(diff(legP,x),xi);
        end
        
        function P = JacobiP(obj, x,alpha,beta,N)
            % function [P] = JacobiP(x,alpha,beta,N)
            % Purpose: Evaluate Jacobi Polynomial of type (alpha,beta) > -1
            %          (alpha+beta <> -1) at points x for order N and returns P[1:length(xp))]
            % Note   : They are normalized to be orthonormal.
            
            % Turn points into row if needed.
            xp = x; dims = size(xp);
            if (dims(2)==1); xp = xp'; end;
            PL = zeros(N+1,length(xp));
            
            % Initial values P_0(x) and P_1(x)
            gamma0 = 2^(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)*...
                gamma(beta+1)/gamma(alpha+beta+1);
            PL(1,:) = 1.0/sqrt(gamma0);
            if (N==0); P=PL'; return; end;
            gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
            PL(2,:) = ((alpha+beta+2)*xp/2 + (alpha-beta)/2)/sqrt(gamma1);
            if (N==1); P=PL(N+1,:)'; return; end;
            
            % Repeat value in recurrence.
            aold = 2/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3));
            
            % Forward recurrence using the symmetry of the recurrence.
            for i=1:N-1
                h1 = 2*i+alpha+beta;
                anew = 2/(h1+2)*sqrt( (i+1)*(i+1+alpha+beta)*(i+1+alpha)*...
                    (i+1+beta)/(h1+1)/(h1+3));
                bnew = - (alpha^2-beta^2)/h1/(h1+2);
                PL(i+2,:) = 1/anew*( -aold*PL(i,:) + (xp-bnew).*PL(i+1,:));
                aold =anew;
            end;
            
            P = PL(N+1,:)';
        end
        
        function dP = dJacobiP(obj, r, alpha, beta, N)
            % function [dP] = GradJacobiP(r, alpha, beta, N);
            % Purpose: Evaluate the derivative of the Jacobi polynomial of type (alpha,beta)>-1,
            %	       at points r for order N and returns dP[1:length(r))]
            
            dP = zeros(length(r), 1);
            if(N == 0)
                dP(:,:) = 0.0;
            else
                dP = sqrt(N*(N+alpha+beta+1))*obj.JacobiP(r(:),alpha+1,beta+1, N-1);
            end;
        end
        
        function [x] = JacobiGL(obj, alpha,beta,N)
            % function [x] = JacobiGL(alpha,beta,N)
            % Purpose: Compute the N'th order Gauss Lobatto quadrature
            %          points, x, associated with the Jacobi polynomial,
            %          of type (alpha,beta) > -1 ( <> -0.5).
            x = zeros(N+1,1);
            if (N==1); x(1)=-1.0; x(2)=1.0; return; end;
            
            [xint,w] = JacobiGQ(alpha+1,beta+1,N-2);
            x = [-1, xint', 1]';
        end
        
        function [x,w] = JacobiGQ(obj, alpha,beta,N)
            
            % function [x,w] = JacobiGQ(alpha,beta,N)
            % Purpose: Compute the N'th order Gauss quadrature points, x,
            %          and weights, w, associated with the Jacobi
            %          polynomial, of type (alpha,beta) > -1 ( <> -0.5).
            
            if (N==0); x(1)= -(alpha-beta)/(alpha+beta+2); w(1) = 2; return; end;
            
            % Form symmetric matrix from recurrence.
            J = zeros(N+1);
            h1 = 2*(0:N)+alpha+beta;
            J = diag(-1/2*(alpha^2-beta^2)./(h1+2)./h1) + ...
                diag(2./(h1(1:N)+2).*sqrt((1:N).*((1:N)+alpha+beta).*...
                ((1:N)+alpha).*((1:N)+beta)./(h1(1:N)+1)./(h1(1:N)+3)),1);
            if (alpha+beta<10*eps); J(1,1)=0.0; end;
            J = J + J';
            
            %Compute quadrature by eigenvalue solve
            [V,D] = eig(J); x = diag(D);
            w = (V(1,:)').^2*2^(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)*...
                gamma(beta+1)/gamma(alpha+beta+1);
        end
        
        function h = DGnflux(obj, f,df,u,strategy)
            %**************************************************************
            % General flux subroutine
            % We use 1 of 4 flux strategies/algorithms optimized for matlab.
            % coded by Manuel Diaz, 2012.12.21 (the end of the world)
            %**************************************************************
            
            % Parameters
            flx = f(u);     % flux value at every point of the domain
            dflx = df(u);   % flux slope at every point of the domain
            
            switch strategy
                case{1} % Roe Flux
                    
                    u_ave = (u(1,:) + u(2,:))/2;    % u @ cells boundaries
                    bool_p = df(u_ave) > 0; % if positive
                    bool_n = df(u_ave) <= 0;  % if negative
                    h = bool_p.*flx(1,:) + bool_n.*flx(2,:);
                    
                case{2} % LF
                    
                    alpha = max(max(abs(dflx)));
                    h = 0.5*(flx(1,:) + flx(2,:) - alpha*(u(2,:)- u(1,:)));
                    
                case{3} % LLF
                    
                    %fluxmat = [dflx(1:nx-1);dflx(2:nx)];
                    beta = max(abs(dflx));
                    h = 0.5*(flx(1,:) + flx(2,:) - beta.*(u(2,:)- u(1,:)));
                    
                case{4} % Upwind Flux
                    
                    % For dflux is constant along the domain!
                    a_p = max(max(dflx - abs(dflx))/2) == [0]; %#ok<NBRAK> % boolen operator for a>0
                    a_n = max(max(dflx + abs(dflx))/2) == [0]; %#ok<NBRAK> % boolen operator for a<0
                    h = a_p*flx(1,:) + a_n*flx(2,:);
                    
                otherwise
                    error('strategy not suported')
            end
        end
        
    end
    
    
end
