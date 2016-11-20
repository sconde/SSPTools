classdef Euler1D < handle
    
    properties
        name = 'Euler-1D';
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
        x; xx;
        y0;
        isSystem = true;
        systemSize = 3;
        domain;
        gamma;
        CFL;
        ProblemType;
        %maxvel;
        BCr;
        BCl;
        BCLvalue;
        BCRvalue;
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
        function obj = Euler1D(varargin)
            % Euler equations of compressible gas dynamics
            p = inputParser;
            p.KeepUnmatched = true;
            
            addParameter(p, 'Gamma', 1.4);
            addParameter(p, 'ProblemType', []);
            addParameter(p, 'N', 100);
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
            
            if ~isempty(p.Results.ProblemType)
                obj.ProblemType = p.Results.ProblemType;
            end
            
            obj.initialize();
        end % Euler1D constructor
        
        function set.xx(obj,fin)
            obj.xx = fin;            
        end
        
        function [p, k, maxvel] = closureModel(obj, density, momentum, energy)
            % ideal gas law
%             density   = u(1:obj.N);
%             momentum  = u(obj.N+1:2*(obj.N));
%             energy   = u(2*(obj.N)+1:3*(obj.N));
            
            p = (obj.gamma - 1)*(energy - 0.5*momentum.^2./density);
            
            c = sqrt(obj.gamma*p./density);
            maxvel = max(c + abs(momentum./density));
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
            
            density   = y(1:nn);
            momentum  = y(nn+1:2*nn);
            energy   = y(2*(nn)+1:3*(nn));
            
            pressure = obj.closureModel(density, momentum, energy);
            
            F = [momentum;                                  % rho*u
                (momentum.^2./density + pressure);          % rho*u^2 + p
                (energy + pressure).*momentum./density];    % (E + p)u
            
        end % flux
        
    end
    
    methods (Access = private )
        
        function initialize(obj)
            
            [r, ru, E] = deal(zeros(obj.N, 1));
            
            if strcmpi(obj.ProblemType, 'sod')
                obj.x = linspace(obj.domain(1), obj.domain(2), obj.N)';
                obj.dx = obj.x(2) - obj.x(1);
                %obj.x = [obj.domain(1):obj.dx:obj.domain(2)]';
                indX = obj.x < 0.5;
                
                r( indX) = 1;
                r(~indX) = 0.125;
                E( indX) = 1/(obj.gamma - 1);
                E(~indX) = 0.1/(obj.gamma - 1);
                
                obj.BCLvalue = [1.0; 0; 2.5];
                obj.BCRvalue = [0.125; 0; 0];
                
                obj.y0 = [r; ru; E];
            elseif strcmpi(obj.ProblemType, 'shock')
                obj.x = linspace(obj.domain(1), obj.domain(2), obj.N)';
                obj.x = linspace(-5, 5, obj.N);
                
                %indX = obj.x < -4;
                for i = 1:obj.N
                    if obj.x(i) < -4
                        rh = 2.857143;
                        u = 2.2629369;
                        p = 10.333333;
                    else
                        rh = 1 + 0.2*sin(pi*obj.x(i));
                        u = 0;
                        p = 1;
                    end
                    r(i) = rh;
                    ru(i) = rh*u;
                    E(i) = p/(obj.gamma - 1) + 0.5*rh*u^2;
                end
                obj.y0 = [r; ru; E];
            end
            
        end % initialize
        
        function [xe, ue] = extend(obj, x, u, h, m, BCl, ul, BCr, ur)
            % extend(x, Q(:,1), h, 1, 'D', 1.0, 'D', 0.125);
            % extend.m : Routine to impose boundary conditions on scalar function by extension
            % Purpose : Extend dependent and independent vectors (x, u), by m cells subject to
            % appropriate boundary conditions
            % BC = 'D' - Dirichlet
            % BC = 'N' - Neumann
            % BC = 'P' - Periodic
            
            xl = min(obj.x);
            xr = max(obj.x);
            N = length(u);
            h = obj.dx;
            m = 1;
            BCl = 'D';
            ul = 1;         % left dirichlet boundary condition
            BCr = 'D';
            ur = 0.125;     % right dirichlet boundary condition
            [xe, ue] = deal(zeros(N + 2*m, 1));
            q = [1:m]';
            
            % Extend x
            xe([m-q+1 N+m+q]) = [xl-q*h xr+q*h];
            xe((m+1):(N+m)) = x(1:N);
            
            % Periodic extension of u
            if (BCl == 'P') | (BCr == 'P')
                ue(m-q+1) = u(N-q);
                ue(N+m+q) = u(q+1);
                ue((m+1):(N+m)) = u(1:N);
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
                ue(N+m+q) = -u(N-q) + 2*ur;
            else
                ue(N+m+q) = u(N-q);
            end
            ue((m+1):(N+m)) = u(1:N);
        end % extend
    end
    
end
