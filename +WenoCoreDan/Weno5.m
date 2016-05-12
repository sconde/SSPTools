classdef Weno5 < WenoCoreDan.Weno
    
    properties
        ep = 1e-8;
        Ngp;
    end
    
    properties (Access = private)
        % Coefficients for the interpolant of each stencil
        stencil_coeffs = [ 2./6., -7./6., 11./6.;
            -1./6.,  5./6., 2./6.;
            2./6.,  5./6., -1./6. ]';
        
        % Weights representing the ideal linear combination of stencils
        ideal_weights = [ 1./10., 6./10., 3./10. ];
        m;
        stencils;
        nonlinear_values;
        wp;
        ISp;
        local_flux_values;
        n_points;
    end
    
    methods
        
        function obj = Weno5(varargin)
            obj = obj@WenoCoreDan.Weno(varargin{:});
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'isSSP', false);
            
            obj.gp = 4;
            obj.order = 5;
            obj.Ngp = obj.nx + 2*obj.gp;
            obj.n_points = obj.nx + 2;
            obj.name = sprintf('WENO%d', obj.order);
            obj.Wp = zeros(obj.Ngp,obj.n_points);
            obj.Wm = zeros(obj.Ngp,obj.n_points);
            obj.WP = zeros(obj.Ngp,obj.n_points);
            obj.WM = zeros(obj.Ngp,obj.n_points);
            obj.m = 3;
            obj.nonlinear_values = zeros(obj.m, obj.m);
            obj.stencils = zeros(obj.m,obj.m);
            obj.wp = zeros(obj.Ngp, obj.n_points);
            
            % Initialize Smoothness Measurements
            obj.ISp = zeros(1,obj.m);
            obj.local_flux_values = zeros(obj.m+(obj.m-1));
        end
        
    end
    
    
    methods
        function [y] = L(obj, t, u)

            % Append ghost points
            u_gp = [ u(end-obj.gp:end-1); u; u(2:obj.gp+1) ];
            
            [fp, fm] = obj.fluxSplit(u_gp, t);
            
            %the differential operator needs the flux?
            u_x = obj.weno_basic(fp, fm)';
            y = u_x;
            
        end
    end
    
    methods (Access = protected)
        
        function [wp] = wenokernel(obj, uLL, uL, u, uR, uRR)
            % this is like the weno function
            
          betaL = (13/12)*(uLL - 2*uL + u).^2 + (1/4)*(uLL - 4*uL + 3*u).^2;
          beta  = (13/12)*(uL - 2*u + uR).^2 + (1/4)*(uL - uR).^2;
          betaR = (13/12)*(u - 2*uR + uRR).^2 + (1/4)*(3*u - 4*uR + uRR).^2;
          
          % do the regularation here
          wL = 0.1./(obj.ep + betaL);
          w  = 0.6./(obj.ep + beta);
          wR = 0.3./(obj.ep + betaR);
          ws = (wL + w + wR);
          wL = wL/ws; 
          wR = wR/ws;
          w = 1 - wL - wR;
          
          % do the reconstruction
          pL = (2*uLL - 7*uL + 11*u)/6;
          p  = (-1*uL + 5*u + 2*uR)/6;
          pR = (2*u + 5*uR - uRR)/6;
          
          wp = wL*pL + w*p + wR*pR;
        end
    end
    
    
end
