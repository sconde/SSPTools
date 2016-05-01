classdef Weno5 < WenoCoreDan.Weno
    
    properties
        
    end
    
    properties (Access = private)
        Ngp;
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
        n_points;
        ISp;
        local_flux_values;
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
            em = obj.em(u_gp);
            
            [fp, fm] = obj.fluxSplit(u_gp, t);
            %the differential operator needs the flux?
            u_x = obj.weno_basic(u_gp, fp, fm)';
            y = u_x;
            
        end
    end
    
    methods (Access = protected)
        
        function [wp] = wenokernel(obj, u_values, u_idx)
            
            persistent stencil_mapping;
            %persistent local_flux_values;
            u_values = u_values';

            
            obj.stencils(:,1) = (1:obj.m)';
            for i=2:obj.m
                obj.stencils(:,i) = obj.stencils(:,i-1) + 1;
            end


            for i=1:obj.m
                stencil_mapping(:, i) = ...
                    sub2ind(size(obj.local_flux_values), i*ones(1,obj.m),...
                    (1:obj.m)+(i-1));
            end
            
            j_absidx = 1;
            
            for j=u_idx
                
                % Get the global indices of the current point.
                window = (j-2:j+2)';
                
                % Get the corresponding function values.
                u = u_values(window);
                
                obj.ISp(1) = 13.0/12.0*(u(1)-2.0*u(2)+u(3))^2 + 0.25*(u(1)-4.0*u(2)+3.0*u(3))^2;
                obj.ISp(2) = 13.0/12.0*(u(2)-2.0*u(3)+u(4))^2 + 0.25*(u(2)-u(4))^2;
                obj.ISp(3) = 13.0/12.0*(u(3)-2.0*u(4)+u(5))^2 + 0.25*(3.0*u(3)-4.0*u(4)+u(5))^2;
                
                % Calculate the normalized weights of each stencil.
                obj.ISp = obj.ideal_weights ./ (obj.ISp + obj.epsilon).^obj.p;
                obj.ISp = obj.ISp ./ norm(obj.ISp,1);
                
                % Apply those weights to each stencil.
                interp_coeffs = bsxfun(@times, obj.stencil_coeffs, obj.ISp);
                
                % Map global indices into our local stencils
                window_indices = window(obj.stencils);
                
                % Update the global flux-interpolation matrix with our local stencils
                obj.local_flux_values(stencil_mapping) = interp_coeffs;
                
                % Update the global flux-interpolation matrix
                obj.wp(window_indices(1):window_indices(end), j_absidx) = ...
                    sum(obj.local_flux_values);
                
                % Move on to the next row of the global flux-interpolation matrix.
                j_absidx = j_absidx + 1;
            end
            
            % Return a correctly-transposed matrix.
            wp = obj.wp';
        end
    end
    
    
end
