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
            
            obj.md = 4; % number of ghost points
            obj.nstart = obj.md + 1;
            obj.remove = obj.md;
            obj.order = 5;
            obj.Ngp = obj.nx + 2*obj.md;
            obj.np = obj.nx + obj.md;
            obj.mpts_inx = (obj.md+1):(obj.nx + obj.md); % ind of good points
            obj.lgpts_inx = (1:obj.md);
            obj.rgpts_inx = (obj.nx+1+obj.md):(obj.nx+2*obj.md);
            obj.xx_inx = [obj.lgpts_inx obj.mpts_inx obj.rgpts_inx];
            obj.dx = 2/(obj.nx); % from [-1 1]
            
            % setup the grid points
            xx = (obj.xx_inx - (obj.md+1))*obj.dx;
            obj.xx = (xx - 1)';
            obj.x = obj.xx(obj.mpts_inx);
            
            % what are these??
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
        function [u_x] = L(obj, t, u)
            
            u = obj.makeu(u);
            % check that we are working with the right size first
            assert(isequal(size(obj.xx), size(u)));
            
            u = obj.applyBC(u);
                       
            em = max(obj.em(u));
            f = obj.f(t, u);
            
            % get the two weno matrices
            [Cn, Dn] = obj.wenokernel(f, u, em);
            ff = f(obj.mpts_inx); uu = u(obj.mpts_inx);
            u_x = 0.5*( Cn* ( ff + em*uu ) + Dn*( ff - em*uu ) )/obj.dx;
        end
    end
    
    methods (Access = protected)
        
        function [C, D] = wenokernel(obj, f, u, em)
            % C <- A == Positive
            % D <- B == Negative
            
            epsilon=10^(-13);
            %global np md nstart remove
                                    
            for im=1:obj.np+obj.md
                fp(im) = 0.5*(f(im) + em*u(im));
                fm(im) = 0.5*(f(im) - em*u(im));
            end            
            
            for im=obj.nstart-obj.md:obj.np+obj.md-1
                dfp(im)= (f(im+1)-f(im) + em*(u(im+1) - u(im)))/2.0;
                dfm(im)= (f(im+1)-f(im) - em*(u(im+1) - u(im)))/2.0;
            end
                        
            for im = obj.nstart-1:obj.np+1
                hh(1,1) = dfp(im-2);
                hh(2,1) = dfp(im-1);
                hh(3,1) = dfp(im);
                hh(4,1) = dfp(im+1);
                hh(1,2) = - dfm(im+2);
                hh(2,2) = - dfm(im+1);
                hh(3,2) = - dfm(im);
                hh(4,2) = - dfm(im-1);
                hh1p(im) = hh(1,1);
                hh2p(im) = hh(2,1);
                hh3p(im) = hh(3,1);
                hh4p(im) = hh(4,1);
                
                hh1m(im) = hh(1,2);
                hh2m(im) = hh(2,2);
                hh3m(im) = hh(3,2);
                hh4m(im) = hh(4,2);
                
                for m1=1:2
                    t1 = hh(1,m1)-hh(2,m1);
                    t2 = hh(2,m1)-hh(3,m1);
                    t3 = hh(3,m1)-hh(4,m1);
                    tt1=(13.0*t1^2 + 3.0*(  hh(1,m1) - 3.0*hh(2,m1))^2);
                    tt2=(13.0*t2^2 + 3.0*(  hh(2,m1) +   hh(3,m1))^2);
                    tt3=(13.0*t3^2 + 3.0*(3.0*hh(3,m1) -   hh(4,m1))^2);
                    tt1=(epsilon+tt1)^2;
                    tt2=(epsilon+tt2)^2;
                    tt3=(epsilon+tt3)^2;
                    s1 = tt2*tt3;
                    s2 = 6.0*tt1*tt3;
                    s3 = 3.0*tt1*tt2;
                    t0 = 1./(s1+s2+s3);
                    s1 = s1*t0;
                    s2 = s2*t0;
                    s3 = s3*t0;
                    w0(m1)=s1;
                    w1(m1)=s2;
                    w2(m1)=s3;
                    if m1 ==1
                        T1p(im) = t1;
                        T2p(im) = t2;
                        T3p(im) = t3;
                        TT1p(im) = tt1;
                        TT2p(im) = tt2;
                        TT3p(im) = tt3;
                        S1p(im) = s1;
                        S2p(im) = s2;
                        S3p(im) = s3;
                        w0p(im) = w0(m1);
                        w1p(im) = w1(m1);
                        w2p(im) = w2(m1);
                    else
                        T1m(im) = t1;
                        T2m(im) = t2;
                        T3m(im) = t3;
                        TT1m(im) = tt1;
                        TT2m(im) = tt2;
                        TT3m(im) = tt3;
                        S1m(im) = s1;
                        S2m(im) = s2;
                        S3m(im) = s3;
                        w0m(im) = w0(m1);
                        w1m(im) = w1(m1);
                        w2m(im) = w2(m1);
                    end
                end
                cp1(im)= (w0(1))/3.0;
                cp2(im)= (-7.0*w0(1) - w1(1) )/6.0;
                cp3(im)= (11.0*w0(1) + 5.0*w1(1) + 2.* w2(1))/6.0;
                cp4(im)= (2.0*w1(1) + 5.0*w2(1) )/6.0;
                cp5(im)= (-w2(1))/6.0;
                cm2(im)= (-w2(2))/6.0;
                cm3(im)= (5.0*w2(2) +2.0*w1(2))/6.0;
                cm4(im)= (2.0*w2(2) + 5.0*w1(2) + 11.0*w0(2))/6.0;
                cm5(im)= (-w1(2) -7.0*w0(2))/6.0;
                cm6(im)= w0(2)/3.0;
            end
                        
            for im=obj.nstart:obj.np
                qp1(im)= cp1(im-1);
                qp2(im)= cp2(im-1) - cp1(im);
                qp3(im)= cp3(im-1) - cp2(im);
                qp4(im)= cp4(im-1) - cp3(im);
                qp5(im)= cp5(im-1) - cp4(im);
                qp6(im) =          - cp5(im);
                qm2(im)= cm2(im-1) ;
                qm3(im)= cm3(im-1) - cm2(im);
                qm4(im)= cm4(im-1) - cm3(im);
                qm5(im)= cm5(im-1) - cm4(im);
                qm6(im)= cm6(im-1) - cm5(im);
                qm7(im)=           - cm6(im);
            end
            
            
            for ii = obj.nstart-obj.remove:obj.np-obj.remove;
                ff(ii)=f(ii+obj.remove);
                uu(ii)=u(ii+obj.remove);
            end
            
            
            % build the matrix:
            for ii = obj.nstart:obj.np;
                for ij=obj.nstart:obj.np;
                    A(ii,ij)= 0;
                    B(ii,ij)= 0;
                end
            end
            
            %work on matrix to make it periodic bc
            for ii=obj.nstart:obj.np;
                A(ii,ii)= qp4(ii);
                B(ii,ii)= qm4(ii);
            end
            
            % sub diagonal
            for ii=obj.nstart+1:obj.np;
                A(ii,ii-1)= qp3(ii);
                B(ii,ii-1)= qm3(ii);
            end;
            
            A(obj.nstart,obj.np) = qp3(obj.nstart);
            B(obj.nstart,obj.np) = qm3(obj.nstart);
            % 2nd sub diagonal
            for ii=obj.nstart+2:obj.np;
                A(ii,ii-2)= qp2(ii);
                B(ii,ii-2)= qm2(ii);
            end;
            
            A(obj.nstart+1,obj.np) = qp2(obj.nstart+1);
            A(obj.nstart,obj.np-1) = qp2(obj.nstart);
            B(obj.nstart+1,obj.np) = qm2(obj.nstart+1);
            B(obj.nstart,obj.np-1) = qm2(obj.nstart);
            
            % 3rd sub diagonal
            for ii=obj.nstart+3:obj.np;
                A(ii,ii-3)= qp1(ii);
                B(ii,ii-3)= 0;
            end;
            
            A(obj.nstart+2,obj.np)   = qp1(obj.nstart+2);
            A(obj.nstart+1,obj.np-1) = qp1(obj.nstart+1);
            A(obj.nstart,obj.np-2)   = qp1(obj.nstart);
            B(obj.nstart+2,obj.np)   = 0;
            B(obj.nstart+1,obj.np-1) = 0;
            B(obj.nstart,obj.np-2)   = 0;
            
            % super diagonal
            for ii=obj.nstart:obj.np-1;
                A(ii,ii+1)= qp5(ii);
                B(ii,ii+1)= qm5(ii);
            end;
            A(obj.np,obj.nstart) = qp5(obj.np);
            B(obj.np,obj.nstart) = qm5(obj.np);
            %
            for ii=obj.nstart:obj.np-2;
                A(ii,ii+2)= qp6(ii);
                B(ii,ii+2)= qm6(ii);
            end;
            
            A(obj.np,obj.nstart+1) = qp6(obj.np);
            A(obj.np-1,obj.nstart) = qp6(obj.np-1);
            B(obj.np,obj.nstart+1) = qm6(obj.np);
            B(obj.np-1,obj.nstart) = qm6(obj.np-1);
            %
            for ii=obj.nstart:obj.np-3;
                B(ii,ii+3)= qm7(ii);
            end;
            B(obj.np,obj.nstart+2)   = qm7(obj.np);
            B(obj.np-1,obj.nstart+1) = qm7(obj.np-1);
            B(obj.np-2,obj.nstart)   = qm7(obj.np-2);
            for ii=obj.nstart:obj.np
                for ij=obj.nstart:obj.np
                    C(ii-obj.remove,ij-obj.remove)=A(ii,ij);
                    D(ii-obj.remove,ij-obj.remove)=B(ii,ij);
                end
            end
        end
    end
    
end
