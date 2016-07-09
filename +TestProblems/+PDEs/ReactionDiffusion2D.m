classdef ReactionDiffusion2D < handle
    
    properties
        name = 'Reaction-Diffusion 2D';
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
        systemSize = 2;
        % parameter for the reaction-diffusion
        a;
        b;
        c;
        d;
        forcingTerm;
    end
    
    properties( Access = private)
        %         ProblemType;
        %         %maxvel;
        %         BCr;
        %         BCl;
        %         BCLvalue;
        %         BCRvalue;
        flagin;
        F;
        X_;
        Y_;
        nn;
        B;
    end
    
    methods
        function obj = ReactionDiffusion2D(varargin)
            % 2D Reaction-Diffusion Systems
            p = inputParser;
            p.KeepUnmatched = true;
            
            addParameter(p, 'a', 1);
            addParameter(p, 'b', 0.5);
            addParameter(p, 'c', 1);
            addParameter(p, 'd', -11);
            addParameter(p, 'SourceTerm',[]);
            
            
            p.parse(varargin{:});
            
            obj.a = p.Results.a;
            obj.b = p.Results.b;
            obj.c = p.Results.c;
            obj.d = p.Results.d;
            obj.forcingTerm = p.Results.SourceTerm;
            
            assert(~isempty(obj.forcingTerm),'Need a forcing Function');
            assert(isequal(length(obj.forcingTerm), obj.systemSize),...
                'num(forcing Function) != size(system)');
                        
        end % ReactionDiffusion2D constructor
        
        function F = flux(obj, t, u)
            x_ = obj.X_(obj.flagin);
            y_ = obj.Y_(obj.flagin);
            
            n_ = length(u); nhalf_ = n_/2;
            u1_ = u(1:nhalf_); u2_ = u(nhalf_+1:end);
            u0 = [u1_(obj.flagin); u2_(obj.flagin)];
            
            f1 = @(t) obj.forcingTerm{1}(t, x_, y_, obj.a, obj.b, obj.c, obj.d);
            f2 =  @(t) obj.forcingTerm{2}(t, x_, y_, obj.a, obj.b, obj.c, obj.d);
            
            u1_ = u0(1:obj.nn);
            u2_ = u0(obj.nn+1:end);
            
            recTerm = [u1.*u2; u1.*u2];
            forcingTerm_ = [f1(t); f2(t)];
            
            F = obj.B*u0 + obj.F.*recTerm + forcingTerm_;
            
        end % flux
        
    end
    
    methods %(Access = private )
        
        function dfdx = initialize(obj, dfdx_obj)
            
            dfdx_obj.D(~dfdx_obj.flagin,:) = [];
            dfdx_obj.D(:, ~dfdx_obj.flagin) = [];
            B_ = blkdiag(obj.a*dfdx_obj.D, obj.b*dfdx_obj.D);
            
            obj.nn = sum(dfdx_obj.flagin);
            e = ones(obj.nn, 1);
            dfdx = dfdx_obj;
            dfdx.D = B_;
            obj.F = [obj.c*e; obj.d*e];
                        
            obj.flagin = dfdx_obj.flagin;
            
            obj.X_ = dfdx.x;
            obj.Y_ = dfdx.y;
            obj.B = B_;
            
            %dfdx.L = @(t, y) obj.flux( t, u);
            keyboard
        end % initialize
        
       
    end
    
end
