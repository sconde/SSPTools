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
        direction;
        flagin;
        f;
        em;
    end
    
    properties (Access = private)
       domainStencil; 
       orderAccuracy;
       stencilSize;
       isInitialized = false;
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
            addParameter(p, 'Direction', 'FD');
            addParameter(p, 'OrderAccuracy', 1);
            addParameter(p, 'Dimension', 1);
            p.parse(varargin{:});
            
            if ~isempty(p.Results.Problem)
                obj.problem = p.Results.Problem;
            end
            
            obj.isSSP = p.Results.isSSP;
            obj.bc = p.Results.bc;
            obj.derivativeOrder = p.Results.derivativeOrder;
            obj.direction = p.Results.Direction;
            obj.orderAccuracy = p.Results.OrderAccuracy;
            obj.stencilSize = obj.orderAccuracy + 1;
            obj.dimN = p.Results.Dimension;
            
            if strcmp(obj.direction,'CD') % centered finite difference
                obj.domainStencil = -obj.orderAccuracy/2:obj.orderAccuracy/2;
            elseif strcmp(obj.direction,'FD') % forward finite difference
                obj.domainStencil = 0:obj.stencilSize-1;
            elseif strcmp(obj.direction,'BD') % backward finite difference
                obj.domainStencil = -(obj.stencilSize-1:-1:0);
            end
            
            
            obj.domain = p.Results.domain;
            obj.nx = p.Results.N;
            
            if ~isempty(obj.domain)
                obj.x = linspace(obj.domain(1),obj.domain(2),obj.nx);
                obj.dx = obj.x(2) - obj.x(1);
            end
            
            
            if ~obj.problem.isSystem
                obj.systemSize = 1;
            else
                
                if ~isa(obj.problem, 'TestProblems.PDEs.ReactionDiffusion2D')
                    obj.dx      = obj.problem.dx;
                    obj.x       = obj.problem.x;
                    obj.nx       = obj.problem.N;
                    obj.domain  = obj.problem.domain;
                    obj.systemSize = obj.problem.systemSize;
                else
                    
                    % initialize the discretizer
                    obj.initialize();
                    
                    % now initialize the 2D systems
                    
                    %call the problem initialize method
                    obj.problem.initialize(obj);
                    %keyboard
                end
            end
            
            obj.initialize();
        end % end constructor
        
    end
    
    methods
        
        function [y] = L(obj, t, y)
            
            % why is this being called??
                        
            yx = zeros(size(y));
            IDX = obj.nx*[0:obj.systemSize-1 ; 1:obj.systemSize]';
            IDX = IDX + [ones(size(IDX(:,1))) zeros(size(IDX(:,1)))];
            fx = obj.problem.f(t,y);
            for i = 1:obj.systemSize
                idx = IDX(i,1); idx_ub = IDX(i,2);
                yx(idx:idx_ub) = obj.D*fx(idx:idx_ub);
            end
            y = yx;
        end % L
    end
    
    methods (Access = private)
        
        function initialize(obj)
            
            if ~obj.isInitialized
                
                obj.nx = numel(obj.x);
                obj.x = obj.x(:);
                
                periodic_ = ~strcmpi('non-periodic',obj.bc);

                if obj.derivativeOrder == 2
                    one = ones(obj.nx , 1);
                    D2 = spdiags([one -2*one one], -1:1, obj.nx , obj.nx);
                    D2(1, 1:4) = [2 -5 4 -1];
                    D2(end, end - 3:end) = [-1 4 -5 2];
                    D2= 1/(obj.dx^obj.derivativeOrder)*D2;
                    T_ = D2;
                else
                    T_ = diffMatix(obj, obj.derivativeOrder, obj.orderAccuracy, ...
                        obj.nx, periodic_, obj.domainStencil)';
                    T_ = T_/obj.dx;
                end


                
                % for 2D case
                if obj.dimN > 1
                   % at the moment, just assume square, structured grid
                   obj.domain = repmat(obj.domain,obj.dimN,1);
                   obj.y = obj.x;
                   obj.dy = obj.dx;
                   
                   [X, Y] = meshgrid(obj.x, obj.y);
                   flagin_ = false(size(X));
                   flagin_(2:end-1, 2:end-1) = true;
                   flagin_ = flagin_(:);
                   I_ = speye(obj.nx);
                   T_ = kron(T_, I_) + kron(I_, T_);
                   obj.x = X;
                   obj.y = Y;
                   obj.flagin = flagin_;
                   obj.isInitialized = true;
                end
                obj.D = T_;
            end
            
        end % initialize
        
        function c = fdcoeffF(obj, k, xbar, x)
            
            % Compute coefficients for finite difference approximation for the
            % derivative of order k at xbar based on grid values at points in x.
            %
            % This function returns a row vector c of dimension 1 by n, where n=length(x),
            % containing coefficients to approximate u^{(k)}(xbar),
            % the k'th derivative of u evaluated at xbar,  based on n values
            % of u at x(1), x(2), ... x(n).
            %
            % If U is a column vector containing u(x) at these n points, then
            % c*U will give the approximation to u^{(k)}(xbar).
            %
            % Note for k=0 this can be used to evaluate the interpolating polynomial
            % itself.
            %
            % Requires length(x) > k.
            % Usually the elements x(i) are monotonically increasing
            % and x(1) <= xbar <= x(n), but neither condition is required.
            % The x values need not be equally spaced but must be distinct.
            %
            % This program should give the same results as fdcoeffV.m, but for large
            % values of n is much more stable numerically.
            %
            % Based on the program "weights" in
            %   B. Fornberg, "Calculation of weights in finite difference formulas",
            %   SIAM Review 40 (1998), pp. 685-691.
            %
            % Note: Forberg's algorithm can be used to simultaneously compute the
            % coefficients for derivatives of order 0, 1, ..., m where m <= n-1.
            % This gives a coefficient matrix C(1:n,1:m) whose k'th column gives
            % the coefficients for the k'th derivative.
            %
            % In this version we set m=k and only compute the coefficients for
            % derivatives of order up to order k, and then return only the k'th column
            % of the resulting C matrix (converted to a row vector).
            % This routine is then compatible with fdcoeffV.
            % It can be easily modified to return the whole array if desired.
            %
            % From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)
            
            
            n = length(x);
            if k >= n
                error('*** length(x) must be larger than k')
            end
            
            m = k;   % change to m=n-1 if you want to compute coefficients for all
            % possible derivatives.  Then modify to output all of C.
            c1 = 1;
            c4 = x(1) - xbar;
            C = zeros(n-1,m+1);
            C(1,1) = 1;
            for i=1:n-1
                i1 = i+1;
                mn = min(i,m);
                c2 = 1;
                c5 = c4;
                c4 = x(i1) - xbar;
                for j=0:i-1
                    j1 = j+1;
                    c3 = x(i1) - x(j1);
                    c2 = c2*c3;
                    if j==i-1
                        for s=mn:-1:1
                            s1 = s+1;
                            C(i1,s1) = c1*(s*C(i1-1,s1-1) - c5*C(i1-1,s1))/c2;
                        end
                        C(i1,1) = -c1*c5*C(i1-1,1)/c2;
                    end
                    for s=mn:-1:1
                        s1 = s+1;
                        C(j1,s1) = (c4*C(j1,s1) - s*C(j1,s1-1))/c3;
                    end
                    C(j1,1) = c4*C(j1,1)/c3;
                end
                c1 = c2;
            end
            
            c = C(:,end)';            % last column of c gives desired row vector
            
        end % fdcoeffF
        
        function T_ = diffMatix(obj, derivativeOrder_, orderAccuracy_, N_, periodic_, domainStencil_)
            
            N_ = N_ + ~periodic_;
            Nzero_ = N_ - orderAccuracy_ - 1;
                        
            stcen_ = obj.fdcoeffF(derivativeOrder_, 0, domainStencil_);
            
            rtop_ = [stcen_ zeros(1,Nzero_)];
            ctop_ = circshift(fliplr(rtop_),[1 1]);
            T_ = toeplitz(ctop_,rtop_);
            T_(end,:) = T_(end - ~periodic_,:);
            
        end % diffMatix
        
    end
    
    
end
