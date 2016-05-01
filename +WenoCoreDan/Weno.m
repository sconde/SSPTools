classdef Weno < handle
    
    properties
        name;
        bc;
        nx;
        domain;
        dx;
        order;
        derivativeOrder;
        x;
        em;
        f;
    end
    
    properties (Access = protected)
        isSSP = false;
        weno_fcn;
        gp;
        epsilon;
        p;
        Wp;
        Wm;
        WP;
        WM;
        fm;
        fp;
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
            addParameter(p,'weno_fcn', @WenoCoreDan.weno_basic);
            p.parse(varargin{:});
            
            obj.isSSP = p.Results.isSSP;
            obj.bc = p.Results.bc;
            obj.derivativeOrder = p.Results.derivativeOrder;
            obj.nx = p.Results.N;
            obj.domain = p.Results.domain;
            obj.epsilon = p.Results.epsilon;
            obj.p = p.Results.p;
            
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
        
        function [fp, fm] = fluxSplit(obj, u, t)
            % Split the numerical flux into positive and negative components
            fu = obj.f(u, t);
            u_em = 0.5*obj.em(u);
            fp = 0.5*fu + u_em;
            fm = 0.5*fu - u_em;
        end
        
        function [u_x] = weno_basic(obj,u, fp, fm)
            
%             % Compute the upwind interpolation
            %obj.Wp(obj.gp-1:length(fp)-obj.gp,:) = ...
                obj.Wp = obj.wenokernel(fp, obj.gp-1:length(fp)-obj.gp)';
                
%             % Compute the downwind interpolation
                obj.Wm = obj.wenokernel(fm(end:-1:1), obj.gp-1:length(fp)-obj.gp+1)';
                %TODO: size(obj.Wp) != size(obj.Wm). Why??
            
                keyboard
            for i=(obj.gp+1:length(u)-obj.gp)
                obj.WP(i,:) = -(obj.Wp(i,:) - obj.Wp(i-1,:))/obj.dx;
                obj.WM(i,:) = -(obj.Wm(i,:) - obj.Wm(i-1,:))/obj.dx;
            end
            keyboard
            
            u_x = obj.WP*fp + obj.WM*fm;
            u_x = u_x(obj.gp+1:length(u)-obj.gp);
        end
        
        function [wp, is_values] = wenokernel(obj,u, idx) end
    end
    
    
end
