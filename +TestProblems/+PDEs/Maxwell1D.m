classdef Maxwell1D < handle
    
    properties
        name = 'Maxwell-1D';
        %f; % flux
        haveExact = false;
        eqn;
        CFL_MAX = 1;
        isLinear = false;
        em;
        prbTtype;
        N;
        L;
        dx;
        x;
        y0;
        isSystem = true;
        systemSize = 2;
        domain;
        gamma;
        CFL;
        ProblemType;
        %maxvel;
        BCr;
        BCl;
        BCLvalue;
        BCRvalue;
        epl;
        epr;
        mul;
        mur;
        ep;
        mu;
        xx;
    end
    
    properties( Access = private)
        %         ProblemType;
        %         %maxvel;
        %         BCr;
        %         BCl;
        %         BCLvalue;
        %         BCRvalue;
    end
    
    methods
        function obj = Maxwell1D(varargin)
            % Euler equations of compressible gas dynamics
            p = inputParser;
            p.KeepUnmatched = true;
            
            addParameter(p, 'Gamma', 1.4);
            addParameter(p, 'ProblemType', []);
            addParameter(p, 'N', 100);
            addParameter(p, 'Epl', 1.0);
            addParameter(p, 'Mul', 1.0);
            addParameter(p, 'Epr', 2.25);
            addParameter(p, 'Mur', 1.0);
            addParameter(p, 'L', 1);
            addParameter(p, 'BCr', 'D');
            addParameter(p, 'BCl', 'D');
            addParameter(p, 'CFL', 0.90);
            
            
            p.parse(varargin{:});
            
            obj.domain = [0 1];
            obj.gamma = p.Results.Gamma;
            obj.N     = p.Results.N + 1;
            obj.L     = p.Results.L;
            obj.dx    = obj.L/obj.N;
            obj.BCr   = p.Results.BCr;
            obj.BCl   = p.Results.BCl;
            obj.CFL   = p.Results.CFL;
            obj.epl   = p.Results.Epl;
            obj.epr   = p.Results.Epr;
            obj.mul   = p.Results.Mul;
            obj.mur   = p.Results.Mur;
            
            if ~isempty(p.Results.ProblemType)
                obj.ProblemType = p.Results.ProblemType;
            end
            
            obj.initialize();
        end % Maxwell1D constructor
        
        function set.xx(obj,fin)
            obj.xx = fin;
            [~, ~, obj.ep, obj.mu] = CavityExact(obj, 0, obj.xx);
            
        end
        
        
        function [Ef, Hf, ep, mu] = CavityExact(obj, t, x)
            % Purpose: set the exact solution to EM cavity problem
            
            xL = length(x);
            n1 = sqrt(obj.epl*obj.mul);
            n2 = sqrt(obj.epr*obj.mur);
            [Ef, Hf, ep, mu] = deal(zeros(xL,1));
            ii = sqrt(-1);
            
            %Compute omega to match coefficients - set initial guess to obtain different solutions
            omega = fzero(@(x) (n1*tan(n2*x) + n2*tan(n1*x)),5);
            
            % Set up exact solution
            A1 = complex(n2*cos(omega*n2)/(n1*cos(n1*omega)));
            A2 = exp(-ii*omega*(n1+n2));
            B1 = A1*exp(-ii*2*n1*omega);
            B2 = A2*exp(ii*2*n2*omega);
            
            for j = 1:xL
                if (x(j) <= 0)
                    A = A1;
                    B = B1;
                    n = n1;
                    ep(j) = obj.epl;
                    mu(j) = obj.mul;
                else
                    A = A2;
                    B = B2;
                    n = n2;
                    ep(j) = obj.epr;
                    mu(j) = obj.mur;
                end
                Eh = (A*exp(ii*n*omega*x(j)) - B*exp(-ii*n*omega*x(j)))*exp(-ii*omega*t);
                Hh = n*(A*exp(ii*n*omega*x(j)) + B*exp(-ii*n*omega*x(j)))*exp(-ii*omega*t);
                Ef(j) = real(Eh);
                Hf(j) = real(Hh);
            end
        end
        
        function [p, k, maxvel] = closureModel(obj, u)
            u;
            p = ones(size(u));
            maxvel = max(1./sqrt(obj.ep.*obj.mu));
            k = obj.CFL*obj.dx/maxvel;
        end
        
        function [xe, ue] = applyBC(obj, u, BCl, ul, BCr, ur)
            % extend.m : Routine to impose boundary conditions on scalar function by extension
            % Purpose : Extend dependent and independent vectors (x, u), by m cells subject to
            % appropriate boundary conditions
            % BC = 'D' - Dirichlet
            % BC = 'N' - Neumann
            % BC = 'P' - Periodic
            
            xl = min(obj.x);
            xr = max(obj.x);
            NN = length(u);
            m = 1;
            [xe, ue] = deal(zeros(NN + 2*m, 1));
            q = [1:m]';
            
            % Extend x
            xe([m-q+1 NN+m+q]) = [xl-q*obj.dx xr+q*obj.dx];
            xe((m+1):(NN+m)) = obj.x(1:NN);
            
            % Periodic extension of u
            if (BCl == 'P') || (BCr == 'P')
                ue(m-q+1) = u(NN-q);
                ue(NN+m+q) = u(q+1);
                ue((m+1):(NN+m)) = u(1:NN);
                return;
            end
            
            % Left extension
            if BCl == 'D'
                ue(m-q+1) = -u(q+1) + 2*ul;
            else
                ue(m-q+1) = u(q+1);
            end
            
            % Right extension
            if BCr == 'D'
                ue(NN+m+q) = -u(NN-q) + 2*ur;
            else
                ue(NN+m+q) = u(NN-q);
            end
            ue((m+1):(NN+m)) = u(1:NN);
            
        end
        
        function F = f(obj, t, y)
            nn = size(y,1)/obj.systemSize;
            
            Ef   = y(1:nn);
            Eh  = y(nn+1:2*(nn));
            
            F = [Ef./obj.ep;            
                Eh./obj.mu];     
            
        end % flux
        
    end
    
    methods (Access = private )
        
        function initialize(obj)
                        
            if strcmpi(obj.ProblemType, 'cavity')
                % Purpose: set the exact solution to EM cavity problem
                obj.x = linspace(obj.domain(1), obj.domain(2), obj.N)';
                [Ef, Hf, obj.ep, obj.mu] = CavityExact(obj, 0, obj.x);
                obj.y0 = [Ef; Hf];

                % % PEC boundary conditions by mirrow principle
                % [xe,Ee] = extend(x,EM(:,1),h,1,’D’,0,’D’,0); 
                % [xe,He] = extend(x,EM(:,2),h,1,’N’,0,’N’,0);
                obj.BCLvalue = [0; 0];
                obj.BCRvalue = [0; 0];

            end
            
        end % initialize
        
       
    end
    
end
