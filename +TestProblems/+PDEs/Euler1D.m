classdef Euler1D < handle
    
    properties
        name = 'Euler-1D';
        f; % flux
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
        systemSize = 3;
        domain;
    end
    
    properties( Access = private)
        gamma;
        ProblemType;
        maxvel;
        BCr;
        BCl;
        BCLvalue;
        BCRvalue;
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
            
            
            p.parse(varargin{:});
            
            obj.domain = [0 1];
            obj.gamma = p.Results.Gamma;
            obj.N     = p.Results.N;
            obj.L     = p.Results.L;
            obj.dx    = obj.L/obj.N;
            obj.BCr   = p.Results.BCr;
            obj.BCl   = p.Results.BCl;
            obj.f     = @(t,u) obj.rhs(u);
            
            if ~isempty(p.Results.ProblemType)
                obj.ProblemType = p.Results.ProblemType;
            end
            
            obj.initialize();
        end % Euler1D constructor
        
        function FU = rhs(obj, u) % acting like the L method?
            
            %Q = [r ru E];
            density   = u(1:obj.N);
            momentum  = u(obj.N+1:2*(obj.N));
            energy    = u(2*(obj.N)+1:3*(obj.N));
            
            if strcmpi(obj.ProblemType, 'sod')
                % Extend data and assign boundary conditions
                [~, re] = obj.applyBC(density,  obj.BCLvalue(1), obj.BCRvalue(1));   % rho extended
                [~, me] = obj.applyBC(momentum, obj.BCLvalue(2), obj.BCRvalue(2));   % momentum extended
                [~, Ee] = obj.applyBC(energy,   obj.BCLvalue(3), obj.BCRvalue(3));   % Energy extension
            end
            
            %Q = [re(2:obj.N+1); me(2:obj.N+1); Ee(2:obj.N+1)];
            Q = [re(2:end); me(2:end); Ee(2:end)]; %TODO: is this correct?
            
            FU = obj.flux(Q);
            
        end % rhs
        
    end
    
    methods (Access = private )
        
        function F = flux(obj, y)
            
            density   = y(1:obj.N+1);
            momentum  = y(obj.N+2:2*(obj.N+1));
            energy   = y(2*(obj.N+1)+1:3*(obj.N+1));

            pressure = obj.closureModel(y);
            
            F = [momentum;                                  % rho*u
                (momentum.^2./density + pressure);          % rho*u^2 + p
                (energy + pressure).*momentum./density];    % (E + p)u
            
        end % flux
        
        function p = closureModel(obj, u)
            % ideal gas law
            density   = u(1:obj.N+1);
            momentum  = u(obj.N+2:2*(obj.N+1));
            energy   = u(2*(obj.N+1)+1:3*(obj.N+1));
            
            p = (obj.gamma - 1)*(energy - 0.5*momentum.^2./density);
        end
        
        function [xe, ue] = applyBC(obj, u, ul, ur)
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
            if (obj.BCl == 'P') || (obj.BCr == 'P')
                ue(m-q+1) = u(NN-q);
                ue(NN+m+q) = u(q+1);
                ue((m+1):(NN+m)) = u(1:NN);
                return;
            end
            
            % Left extension
            if obj.BCl == 'D'
                ue(m-q+1) = -u(q+1) + 2*ul;
            else
                ue(m-q+1) = u(q+1);
            end
            
            % Right extension
            if obj.BCr == 'D'
                ue(NN+m+q) = -u(NN-q) + 2*ur;
            else
                ue(NN+m+q) = u(NN-q);
            end
            ue((m+1):(NN+m)) = u(1:NN);
            
        end
        
        function initialize(obj)
            
            [r, ru, E] = deal(zeros(obj.N+1, 1));
            
            if strcmpi(obj.ProblemType, 'sod')
                obj.x = [obj.domain(1):obj.dx:obj.domain(2)]';
                indX = obj.x < 0.5;
                
                r( indX) = 1;
                r(~indX) = 0.125;
                E( indX) = 1/(obj.gamma - 1);
                E(~indX) = 0.1/(obj.gamma - 1);
                
                obj.BCLvalue = [1.0; 0; 2.5];
                obj.BCRvalue = [0.125; 0; 0];
                
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
