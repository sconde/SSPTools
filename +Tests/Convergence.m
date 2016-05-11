classdef Convergence < Tests.Test
    
    properties
        refinement_type = 'time';
        CFL;
        t;
        dudt; % should be the integrator
        verbose;
        refine_level;
        Tfinal;
        referenceSolution;
        epss;
        L2Error;
        L1Error;
        LinfError;
        DT;
        problemName;
    end
    
    properties ( Access = private)
        y0;
        dx;
        u0;
        L2Error_ = [];
        L1Error_ = [];
        LinfError_ = [];
    end
    
    methods
        function obj = Convergence(varargin)  %constructor
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('refinement_type', []);
            p.addParameter('CFL', []);
            p.addParameter('t', 0);
            p.addParameter('integrator', []);
            p.addParameter('verbose', 'none'); %TODO: decide if type should be logical instead
            p.addParameter('refine_level', 5);
            p.addParameter('Tfinal', 0.1);
            p.addParameter('L2Error', true);
            p.addParameter('L1Error', true);
            p.addParameter('LinfError', true);
            p.addParameter('NonlinearEps', 1e-14);
            p.parse(varargin{:});
            
            obj = obj@Tests.Test(varargin{:});
            obj.refinement_type = p.Results.refinement_type;
            obj.refine_level = p.Results.refine_level;
            %             obj.getL1Error = p.Results.L1Error;
            %             obj.getL2Error = p.Results.L2Error;
            %             obj.
            
            if ~isempty(p.Results.CFL)
                obj.CFL = p.Results.CFL;
            else
                obj.CFL = 0.2*(1/2).^(1:obj.refine_level);
            end
            
            obj.t = p.Results.t;
            obj.verbose = p.Results.verbose;
            obj.Tfinal = p.Results.Tfinal;
            obj.refinement_type = p.Results.refinement_type;
            obj.epss = p.Results.NonlinearEps;
            
            %TODO : a better way to test for the problem
            obj.dudt = p.Results.integrator;
            obj.DT = obj.dudt.dfdx.dx*obj.CFL;
            
            
            if isa(obj.dudt,'SSPTools.Steppers.IMEXRK')
                
                %if IMEX: use ODE45 to get the solution vector
                ode45_options = odeset('RelTol',obj.epss,'AbsTol',obj.epss);
                warning('off');
                
                odefunc = @(t,y) obj.dudt.dfdx.L(obj.dudt.ExplicitProblem.f(t,y))...
                    + obj.dudt.dgdx.L(obj.dudt.ImplicitProblem.f(t,y));
                
                [~,sol] = ode45(@(t,y) odefunc(t, y),[0 obj.Tfinal],...
                    obj.dudt.y0(obj.dudt.dfdx.x),ode45_options);
                
                sol = sol(end,:); sol = sol(:);
                obj.referenceSolution = sol;
            elseif isa(obj.dudt.ExplicitProblem,'TestProblems.PDEs.LinearAdvection')
                obj.referenceSolution = ...
                    obj.dudt.y0(obj.dudt.dfdx.x - obj.dudt.ExplicitProblem.a * obj.Tfinal);
                
            elseif ~isa(obj.dudt, 'SSPTools.Steppers.IMEXRK') %is not using IMEXRK
                % use ODE45 to get the solution vector
                
                ode45_options = odeset('RelTol',obj.epss,'AbsTol',obj.epss);
                warning('off');
                odefunc = @(t,y) obj.dudt.dfdx.L(obj.dudt.ExplicitProblem.f(t,y));
                [~,sol] = ode45(@(t,y) odefunc(t, y),[0 obj.Tfinal],...
                    obj.dudt.y0(obj.dudt.dfdx.x),ode45_options);
                sol = sol(end,:); sol = sol(:);
                obj.referenceSolution = sol;
                
            else
                error('not yet implemented');
            end
            
            if isa(obj.dudt,'SSPTools.Steppers.IMEXRK')
                obj.problemName = sprintf('%s (Exp) + %s (Imp)',...
                    obj.dudt.ExplicitProblem.name, ...
                    obj.dudt.ImplicitProblem.name);
            else
                %TODO: is this working?
                obj.problemName = sprintf('%s',...
                    obj.dudt.ExplicitProblem.name);
            end
        end
        
        function run(obj, varargin)
            % just the temporal refinement for now
            
            % this is where I'm taking the steps
            for dt = obj.DT
                t_ = 0;
                obj.dudt.resetInitCondition();
                
                while t_ < obj.Tfinal
                    obj.dudt.takeStep(dt);
                    t_ = t_ + dt;
                    dt = min(dt, obj.Tfinal - t_);
                end
                [t_, y] = obj.dudt.getState();
                assert(isequal(obj.Tfinal, t_),'should be calculating error at tfinal');
                
                err = abs(y - obj.referenceSolution);
                obj.L1Error_ = [obj.L1Error_; norm(err, 1)];
                obj.L2Error_ = [obj.L2Error_; norm(err, 2)];
                obj.LinfError_ = [obj.LinfError_; norm(err, 'inf')];
            end
            obj.L1Error = obj.L1Error_;
            obj.L2Error = obj.L2Error_;
            obj.LinfError = obj.LinfError_;
        end
        
        function name = getStepperName(obj)
            name = obj.dudt.name;
        end
        
        function order = getOrder(obj,Err)
            
            if strcmpi(Err, 'l2')
                err = obj.L2Error;
            elseif strcmpi(Err, 'l1')
                err = obj.L1Error;
            elseif strcmpi(Err, 'linf')
                err = obj.LinfError;
            end
            obsOrder = polyfit(log(obj.DT(:)),log(err),1);
            order = obsOrder(1);
        end
        
        
        function complete(obj)
            % print the result of the test
            fid = 1;
            getOrder = @(err, dt) log(err(2:end)./err(1:end-1))'./log(dt(2:end)./dt(1:end-1));
            
            Order_inf = [nan getOrder(obj.LinfError, obj.DT)]';
            Order_l1 = [nan getOrder(obj.L1Error, obj.DT)]';
            Order_l2 = [nan getOrder(obj.L2Error, obj.DT)]';
            
            fprintf(fid, '%s (Tfinal = %4.3f)\n\n', obj.problemName, obj.Tfinal);
            order_info= [obj.DT' Order_inf Order_l1 Order_l2];
            fprintf(fid, '%8s \t %9s \t %6s \t %6s\n','DT','L-inf','L1','L2');
            fprintf(fid, '%6.2e \t %6.5f \t %6.5f \t %6.5f\n',order_info(:,:).');
        end
        
    end
    
    
    methods %(Access = protected)
        %         function log(obj, varargin) end
    end
    
    methods (Access = private)
        
        %         function refine_time_only(obj, varargin)
        %
        % %         end
        %
        %         function dy = myodefunSpectral(t, y)
        %             t = 0;
        %             dy = lin_func(y) + nonlinear_func(y);
        %         end
    end
    
end
