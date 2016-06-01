classdef Weno < handle
    
    properties
        name;
        bc;
        nx;
        domain;
        dx;
        order;
        derivativeOrder = 1;
        x;
        em;
        f;
        xx; % x with ghost points
        mpts_inx;  % good point
        lgpts_inx; % left ghost point index
        rgpts_inx; % right ghost point index
        xx_inx; % all the indeces
        problem;
    end
    
    properties (Access = protected)
        isSSP = false;
        weno_fcn;
        md;
        epsilon;
        p;
        Wp;
        Wm;
        WP;
        WM;
        fm;
        fp;
        nstart;
        remove;
        np;
    end
    
    methods
        
        function obj = Weno(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'isSSP', false);
            addParameter(p, 'derivativeOrder', 1);
            addParameter(p, 'bc', 'periodic');
            addParameter(p, 'N', 10);
            addParameter(p, 'domain', [-1 1]);
            addParameter(p,'em', []);
            addParameter(p,'epsilon', 1e-16);
            addParameter(p,'p', 2);
            addParameter(p,'kernel', []);
            addParameter(p, 'Problem', []);
            addParameter(p,'weno_fcn', @WenoCoreDan.weno_basic);
            p.parse(varargin{:});
            
            obj.isSSP = p.Results.isSSP;
            obj.bc = p.Results.bc;
            obj.derivativeOrder = p.Results.derivativeOrder;
            obj.nx = p.Results.N;
            obj.domain = p.Results.domain;
            obj.epsilon = p.Results.epsilon;
            obj.p = p.Results.p;
            
            if ~isempty(p.Results.Problem)
                obj.problem = p.Results.Problem;
            end
            
            if ~isempty(obj.domain)
                obj.x = linspace(obj.domain(1),obj.domain(2),obj.nx);
                obj.dx = min(diff(obj.x));
                obj.nx = numel(obj.x);
                obj.x = obj.x(:);
            end
        end
        
        function set.f(obj,fin)
            if isa(fin, 'function_handle')
                obj.f = fin;
            end
        end
        
        function set.em(obj, val)
            if isa(val, 'function_handle')
                obj.em = val;
            end
        end
        
        function [y] = L(obj, u)
            
        end
        
    end
    
    methods (Access = protected)
        
        function u = applyBC(obj, u)
            u(obj.lgpts_inx) = u(obj.mpts_inx(end-obj.md+1:end));
            u(obj.rgpts_inx) = u(obj.mpts_inx(1:obj.md));
        end
        
        function [fp, fm] = fluxSplit(obj, u, t)
            % Split the numerical flux into positive and negative components
            fu = obj.f(u, t);
            u_em = 0.5*obj.em(u);
            fp = 0.5*fu + u_em;
            fm = 0.5*fu - u_em;
            
        end
        
        function unew = makeu(obj, u)
            
            unew = zeros(obj.np+obj.md,1);
            
            for i = obj.nstart:obj.np
                unew(i)= u(i-obj.remove);
            end;
            % periodic boundaries
            for i = 1:obj.md
                unew(obj.np+i) = unew(obj.nstart+i-1);
            end;
            for i = 1:obj.md
                unew(obj.nstart-i)= unew(obj.np+1-i);
            end;
        end
            
            function [u_x] = weno_basic(obj,fp, fm)
                
                % Compute the upwind interpolation
                %obj.Wp(obj.md-1:length(fp)-obj.md,:) = ...
                n = obj.Nmd;
                uInd = obj.md+1:n-obj.md;
                obj.Wp = obj.wenokernel(fp(uInd-2), fp(uInd-1), fp(uInd), fp(uInd+1), fp(uInd+2));
                
                % Compute the downwind interpolation
                obj.Wm = obj.wenokernel(fm(uInd+2), fm(uInd+1), fm(uInd), fm(uInd-1), fm(uInd-2));
                
                error('Not yet finished');
                
                u_x = obj.WP*fp + obj.WM*fm;
                u_x = u_x(obj.md+1:length(u)-obj.md);
            end
            
            function [wp, is_values] = wenokernel(obj,u, idx) end
        end
        
        
    end
