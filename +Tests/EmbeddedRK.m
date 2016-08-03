classdef EmbeddedRK < Tests.Test
    
    properties
        refinement_type = 'eps';
        CFL;
        t;
        dudt; % should be the integrator
        verbose;
        refine_level;
        Tfinal;
        referenceSolution;
        epss;
        NStep;
        NAccept;
        NReject;
        NFunEval;
        DT;
        problemName;
        exactSol;
        TOL; % range of tolerance to test against
        TOL_ = [];
        FNC_ = [];
        STEP_ = [];
        ACCP_ = [];
        REJCT_ = [];
    end
    
    properties ( Access = private)
        y0;
        dx;
        u0;
        fcn_ = [];
        step_ = [];
        accpt_ = [];
        reject_ = [];
        isExactSolSet = false;
        dt_ ;
    end
    
    methods
        function obj = EmbeddedRK(varargin)  %constructor
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('refinement_type', []);
            p.addParameter('t', 0);
            p.addParameter('integrator', []);
            p.addParameter('verbose', 'none'); %TODO: decide if type should be logical instead
            p.addParameter('Tolerance',[]);
            p.addParameter('dt', 0.1);
            p.parse(varargin{:});
            
            obj = obj@Tests.Test(varargin{:});
            obj.refinement_type = p.Results.refinement_type;
            obj.TOL = p.Results.Tolerance;
            obj.name = 'Embedded RK Tests';
            obj.STEP_ = nan(size(obj.TOL));
            obj.FNC_ = nan(size(obj.TOL));
            obj.ACCP_ = nan(size(obj.TOL));
            obj.REJCT_ = nan(size(obj.TOL));
            
            %TODO : a better way to test for the problem
            obj.dudt = p.Results.integrator;    
            
            obj.t = p.Results.t;
            obj.verbose = p.Results.verbose;
            obj.Tfinal = obj.dudt.tFinal;
            obj.refinement_type = p.Results.refinement_type;
            obj.dt_ = p.Results.dt;
                        
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
        
        function run(obj)
            
            for i = 1:numel(obj.TOL)
                tol = obj.TOL(i);
                
                obj.singleRun(tol);
                
                obj.STEP_(i) = obj.dudt.nstep;
                obj.FNC_(i) = obj.dudt.nfcn;
                obj.ACCP_(i) = obj.dudt.naccpt;
                obj.REJCT_(i) = obj.dudt.nrejct;
                
            end
            
            obj.NAccept = obj.ACCP_;
            obj.NReject = obj.REJCT_;
            obj.NFunEval = obj.FNC_;
            obj.NStep = obj.STEP_;
            
        end
        
        function singleRun(obj, tol)
            % just the temporal refinement for now

            T = []; DT = []; Y = []; ERR = []; badDT = [];

			obj.dudt.resetInitCondition(tol);
            [t, y, dt, err, bad_dt] = obj.dudt.getState();
            
            T = [T; t]; DT = [DT; dt]; Y = [Y y]; ERR = [ERR; err]; badDT = [badDT; bad_dt];
                        
            while t < obj.dudt.tFinal
                
                obj.dudt.takeStep(dt);
                [t, y, nextDt, err, bad_dt] = obj.dudt.getState();
                dt = min(nextDt, obj.dudt.tFinal - t);
                T = [T; t]; DT = [DT; dt]; Y = [Y y]; ERR = [ERR; err]; badDT = [badDT; bad_dt];
            end

            %obj.dudt.summary();     

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
        
        function err = getError(obj,Err)
            
            if strcmpi(Err, 'l2')
                err = obj.L2Error;
            elseif strcmpi(Err, 'l1')
                err = obj.L1Error;
            elseif strcmpi(Err, 'linf')
                err = obj.LinfError;
            end
        end
        
        function dt = getDT(obj)
            dt = obj.DT;
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
