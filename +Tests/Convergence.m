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

    end
    
    properties ( Access = private)
        y0;
        DT;
        dx;
        u0;
        L2Error;
        L1Error;
        LinfError;
    end
    
    methods
        function obj = Convergence(varargin)  %constructor
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParamValue('refinement_type', []);
            p.addParamValue('CFL', []);
            p.addParamValue('t', 0);
            p.addParamValue('integrator', []);
            p.addParamValue('verbose', 'none'); %TODO: decide if type should be logical instead
            p.addParamValue('refine_level', 5);
            p.addParamValue('Tfinal', 0.1);
            p.addParamValue('L2Error', true);
            p.addParamValue('L1Error', true);
            p.addParamValue('LinfError', true);
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
            
            %TODO : a better way to test for the problem
            obj.dudt = p.Results.integrator;
            obj.DT = obj.dudt.dfdx.dx*obj.CFL;
            
            odefunc = @(t,y) obj.dudt.dfdx.L(obj.dudt.ExplicitProblem.f(t,y));
                                    
            if isa(obj.dudt.ExplicitProblem,'TestProblems.PDEs.LinearAdvection')
                obj.referenceSolution = obj.dudt.y0(obj.dudt.dfdx.x - obj.dudt.ExplicitProblem.a * obj.Tfinal);
            elseif ~isa(obj.dudt, 'SSPTools.Steppers.IMEXRK')
                % use ODE45 to get the solution vector
                % is not using IMEXRK
                epss = 1e-14;
                ode45_options = odeset('RelTol',epss,'AbsTol',epss);
                [~,sol] = ode45(@(t,y) odefunc(t, y),[0 obj.Tfinal],...
                    obj.dudt.y0(obj.dudt.dfdx.x),ode45_options);
                sol = sol(end,:); sol = sol(:);
                obj.referenceSolution = sol;
                
            else
                error('not yet implemented');
            end
        end
        
        function run(obj, varargin)
            % just the temporal refinement for now
            if strcmpi(obj.refinement_type, 'time')
                obj.refine_time_only(varargin{:});
            end
            
            L1Error = [];
            L2Error = [];
            LinfError = [];
            % this is where I'm taking the steps
            for dt = obj.DT
                t = 0;
                obj.dudt.resetInitCondition();
                
                while t < obj.Tfinal
                    obj.dudt.takeStep(dt);
                    t = t + dt;
                    dt = min(dt, obj.Tfinal - t);
                end
                [t, y] = obj.dudt.getState();
                assert(isequal(obj.Tfinal, t),'should be calculating error at tfinal');
                
%                 clf;
%                 plot(obj.dudt.dfdx.x,y,'-r')
%                 hold on
%                 plot(obj.dudt.dfdx.x,obj.referenceSolution,'-k')
%                 keyboard
                err = abs(y - obj.referenceSolution); %clear('y');
                L1Error = [L1Error; norm(err, 1)];
                L2Error = [L2Error; norm(err, 2)];
                LinfError = [LinfError; norm(err, 'inf')];
            end
            obj.L1Error = L1Error;
            obj.L2Error = L2Error;
            obj.LinfError = LinfError;
        end
        
        function [ output ] = run_test(varargin)
            
        end
        
        
        function complete(obj)
            % print the result of the test
            fid = 1;
            getOrder = @(err, dt) log(err(2:end)./err(1:end-1))'./log(dt(2:end)./dt(1:end-1));
            
            Order_inf = [nan getOrder(obj.LinfError, obj.DT)]';
            Order_l1 = [nan getOrder(obj.L1Error, obj.DT)]';
            Order_l2 = [nan getOrder(obj.L2Error, obj.DT)]';

%TODO: Do I need these?            
%             order_l2 = polyfit(log(obj.DT(:)), log(Err_l2(:)), 1);
%             order_l1 = polyfit(log(obj.DT(:)), log(Err_l1(:)), 1);
%             order_inf = polyfit(log(obj.DT(:)), log(Err_inf(:)), 1);
            
            order_info= [obj.DT' Order_inf Order_l1 Order_l2];
            fprintf(fid, '%8s \t %9s \t %6s \t %6s\n','DT','L-inf','L1','L2');
            fprintf(fid, '%6.2e \t %6.5f \t %6.5f \t %6.5f\n',order_info(:,:).');
            
%             fprintf(fid, '\nPolyfit-Order:\n');
%             fprintf(fid, 'Order:\t L_inf = %5.3f,\t L_1 = %5.3f, \t L_2 = %5.3f', order_inf(1), order_l1(1), order_l2(1));
%             
keyboard
        end
        
%         function order = getOrder(err, dt)
%             order = log(err(2:end)./err(1:end-1))'./log(dt(2:end)./dt(1:end-1));
%         end
        
    end
    
    
    methods %(Access = protected)
        function log(obj, varargin) end
    end
    
    methods (Access = private)
        
%         function order = getOrder(obj,err, dt)
%             order = log(err(2:end)./err(1:end-1))'./log(dt(2:end)./dt(1:end-1));
%         end
        
        function refine_time_only(obj, varargin)
            % Run the test.
            
            results = {};
            completed_problems = {};
            
            obj.log('Convergence Test\n');
            obj.log('Problem: %s\n', obj.problem_template.repr());
            obj.log('Time Stepping Method: %s\n', obj.problem_template.integrator.repr());
            obj.log('Spatial Discretization: %s\n', obj.problem_template.discretizer.repr());
            obj.log('\n');
            
            for i=1:numel(obj.refinements)
                
                n = 100;
                dt = obj.refinements(i);
                
                problem = obj.problem_template.copy();
                problem.setup_problem(n);
                
                obj.log('Testing %i dt=%g\n', n, dt);
                
                problem.approximate(obj.t, 'dt', dt, 'verbose', false);
                l1error = problem.error_norm(1);
                l2error = problem.error_norm(2);
                linferror = problem.error_norm(inf);
                dx = min(diff(problem.x));
                
                results{end+1} = struct('N', n,...
                    'dt', dt,...
                    'dx', dx,...
                    'u', problem.u,...
                    'x', problem.x,...
                    't', problem.t,...
                    'exact', problem.get_exact_solution(), ...
                    'pointwise', problem.calculate_error(),...
                    'l1error', l1error,...
                    'l2error', l2error,...
                    'linferror', linferror);
                
                completed_problems{end+1} = problem;
            end
            
            results = [ results{:} ];
            
            for i=2:numel(results)
                refinement = log10(results(i).dt / results(i-1).dt);
                results(i).l2order  = log10( results(i).l2error / results(i-1).l2error) / refinement;
                results(i).l1order  = log10( results(i).l1error / results(i-1).l1error) / refinement;
                results(i).linforder  = log10( results(i).linferror / results(i-1).linferror) / refinement;
            end
            
            obj.results = results;
            obj.completed_problems = [ completed_problems{:} ];
        end
        
        function dy = myodefunSpectral(t, y)    
            dy = lin_func(y) + nonlinear_func(y);  
        end
    end
    
end
