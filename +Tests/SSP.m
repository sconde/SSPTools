classdef SSP < Tests.Test
    
    properties
        cfl;
        t;
        dudt; % should be the integrator
        verbose;
        Tfinal;
        ssp;
        theortical_r;
        TV;
        log10VV;
        CFL;
    end
    
    properties ( Access = private)
        y0;
        DT;
        dx;
        u0;
        problemName;
        initTV;
        numViolation;
        Steps;
        testTVB;
        testTVD;
        TVB;
        TVD;
        CFLMAX;
        cflRefine;
        sspCoef; % theorical ssp coef
        startingr;
        lambda;
        cfl_refinement;
        acc;
        tol;
    end
    
    
    methods
        function obj = SSP(varargin)  %constructor
            obj = obj@Tests.Test(varargin{:});
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('cfl', 0.2);
            p.addParameter('t', 0);
            p.addParameter('integrator', []);
            p.addParameter('verbose', 'none'); %TODO: decide if type should be logical instead
            p.addParameter('Tfinal', 0.1);
            p.addParameter('NumberAllowedViolation', 10);
            p.addParameter('Steps', 20);
            p.addParameter('TVB', false);
            p.addParameter('TVD', false);
            p.addParameter('CFLMAX', 1.5);
            p.addParameter('CFLRefinement',0.1);
            p.addParameter('Tolerance', 1e-5);
            p.addParameter('r', 1);
            p.parse(varargin{:});
            
            %obj = obj@Tests.Test(varargin{:});
            
            obj.cfl = p.Results.cfl;
            obj.t = p.Results.t;
            obj.verbose = p.Results.verbose;
            obj.Tfinal = p.Results.Tfinal;
            obj.numViolation = p.Results.NumberAllowedViolation;
            obj.Steps = p.Results.Steps;
            obj.testTVB = p.Results.TVB;
            obj.testTVD = p.Results.TVD;
            obj.CFLMAX = p.Results.CFLMAX;
            obj.cflRefine = p.Results.CFLRefinement;
            obj.theortical_r = p.Results.r;
            obj.tol = p.Results.Tolerance;
            
            %TODO: is this correct?
            %assert(any([obj.testTVD obj.testTVB]),'Must choice TVD or TVB');
            
            
            %TODO : a better way to test for the problem
            obj.dudt = p.Results.integrator;
            obj.DT = obj.dudt.dfdx.dx*obj.cfl;
            [~, y] = obj.dudt.getState();
            obj.initTV = obj.calcTV(y);
            
            if obj.dudt.isSSP
                obj.sspCoef = obj.dudt.r;
            end
            
            obj.initialize();
            
        end
        
        function initialize(obj)
            obj.startingr = max(obj.theortical_r/2,.005);
            %obj.sp = min(obj.theortical_r/10,.005);
            obj.lambda = linspace(obj.startingr, obj.theortical_r + 0.5, 10);
            obj.cfl_refinement = max(diff(obj.lambda));
        end
        
        
        function run(obj, varargin)
            
            L = [];
            VV = [];
            V = 0;
            coarse = 1;
            VV_ = [];
            lambda_ = obj.lambda;
                     
            obj.acc = -5;
            numRefinement = 0;
            while obj.cfl_refinement > obj.tol
                Violation_ = obj.runRange(lambda_);
                
                L = [L,lambda_];
                VV_ = [VV_, Violation_];
                
                
                % find the first violation
                ind = find(log10(Violation_) > obj.acc, 1, 'first');
                
                if isempty(ind);  %If the observed CFL is outside original range
                    maxL = max(lambda_);
                    ind = find(lambda_ == maxL, 1, 'first');
                    lambda_ = linspace(maxL - obj.cfl_refinement, maxL + 2*obj.cfl_refinement,10);
                else

                    Ltemp = sort(L);
                    ind_ = find(Ltemp == lambda_(ind),1,'first');
                    
                    if ind_ > 1
                     lambda_ = sort(linspace(Ltemp(ind_ - 1), Ltemp(ind_), 5))
                    else
                        keyboard
                        newL = Ltemp(ind_);
                        lambda_ = linspace(newL-2*obj.cfl_refinement,newL,10); % need to make sure this is all positive
                    end
                end
                obj.cfl_refinement = max(abs(diff(lambda_)));
                obj.acc = min(obj.acc, floor(log10(obj.cfl_refinement)) - 1) ;
                acct = obj.acc;
                numRefinement = numRefinement + 2;
                ref = obj.cfl_refinement;
            end
            
            obj.CFL = L;
            obj.TV = VV_;
            
            %first calculate the SSP
            obj.calculateSSP();
        end
        
        function Violation_ = runRange(obj, lambda)
            tvdFun = @(u) sum([abs(diff(u)); abs((u(1)-u(end)))]);
            
            Violation_ = inf(1, length(lambda));
            nSteps = obj.Steps;
            for i = 1:numel(lambda);
                dt = lambda(i)*obj.dudt.dfdx.dx;
                
                try
                    obj.dudt.resetInitCondition();
                    [~, y_] = obj.dudt.getState();
                    TV_ = nan(1, nSteps);
                    TV_(1) = tvdFun(y_);
                    
                    for tt = 2:nSteps
                        obj.dudt.takeStep(dt);
                        [~, y_] = obj.dudt.getState();
                        TV_(tt) = tvdFun(y_);
                    end
                    maxDiff = max([diff(TV_),1e-15]);
                    Violation_(i) = maxDiff;
                    log10(maxDiff); %TODO: delete this line

                    if log10(maxDiff) > obj.acc
                        break
                    end
                    
                catch err
                    keyboard
                    break
                end
            end
        end
                    
        function [ output ] = run_test(varargin) end
        
        function plotSolution(obj)
            
            %now plot the solution
            plot(obj.CFL, obj.log10VV, 'kx', 'markersize', 8);
        end
        
        function calculateSSP(obj)
            obj.log10VV = log10(obj.TV);
            extremeInd = obj.log10VV >= 0;
            obj.log10VV(extremeInd) = 0;
            
            % get the index of last stable cfl (this is the observed SSP)
            indGood = obj.log10VV <= -14;
            goodSSPInd = find(indGood, 1, 'last');
            obj.ssp = obj.CFL(goodSSPInd);
        end
    end
    
    methods %(Access = protected)
        function log(obj, varargin) end
        
        function tv = calcTV(obj, u)
            tv = sum([abs([diff(u); (u(end)-u(1))])]);
        end
    end
    
end
