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
    end
    
    properties ( Access = private)
        y0;
        DT;
        dx;
        u0;
        L2Error;
        L1Error;
        LinfError;
        problemName;
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
                %if IMEX
                %TODO: implement reference solution for IMEXRK
                % use ODE45 to get the solution vector
                ode45_options = odeset('RelTol',obj.epss,'AbsTol',obj.epss); warning('off');
                odefunc = @(t,y) obj.dudt.dfdx.L(obj.dudt.ExplicitProblem.f(t,y)) + obj.dudt.dgdx.L(obj.dudt.ImplicitProblem.f(t,y));
                [~,sol] = ode45(@(t,y) odefunc(t, y),[0 obj.Tfinal],...
                    obj.dudt.y0(obj.dudt.dfdx.x),ode45_options);
                sol = sol(end,:); sol = sol(:);
                obj.referenceSolution = sol;
                %error('not yet implemented');
            elseif isa(obj.dudt.ExplicitProblem,'TestProblems.PDEs.LinearAdvection')
                obj.referenceSolution = obj.dudt.y0(obj.dudt.dfdx.x - obj.dudt.ExplicitProblem.a * obj.Tfinal);
            elseif ~isa(obj.dudt, 'SSPTools.Steppers.IMEXRK')
                % use ODE45 to get the solution vector
                % is not using IMEXRK
                ode45_options = odeset('RelTol',obj.epss,'AbsTol',obj.epss); warning('off');
                odefunc = @(t,y) obj.dudt.dfdx.L(obj.dudt.ExplicitProblem.f(t,y));
                [~,sol] = ode45(@(t,y) odefunc(t, y),[0 obj.Tfinal],...
                    obj.dudt.y0(obj.dudt.dfdx.x),ode45_options);
                sol = sol(end,:); sol = sol(:);
                obj.referenceSolution = sol;
                
            else
                error('not yet implemented');
            end
            
            obj.problemName = sprintf('%s (Exp) - %s (Imp)',...
                obj.dudt.ExplicitProblem.name, ...
                obj.dudt.ImplicitProblem.name);
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

            fprintf(fid, '%s (Tfinal = %4.3f)\n\n', obj.problemName, obj.Tfinal);
            order_info= [obj.DT' Order_inf Order_l1 Order_l2];
            fprintf(fid, '%8s \t %9s \t %6s \t %6s\n','DT','L-inf','L1','L2');
            fprintf(fid, '%6.2e \t %6.5f \t %6.5f \t %6.5f\n',order_info(:,:).');
        end
        
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
