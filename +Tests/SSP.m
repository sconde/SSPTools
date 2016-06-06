classdef SSP < Tests.Test
    
    properties
        cfl;
        t;
        dudt; % should be the integrator
        verbose;
        Tfinal;
        ssp;
    end
    
    properties ( Access = private)
        y0;
        DT;
        dx;
        u0;
        problemName;
        initTV;
        TV;
        numViolation;
        Steps;
        testTVB;
        testTVD;
        TVB;
        TVD;
        CFLMAX;
        cflRefine;
        sspCoef; % theorical ssp coef
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
            p.addParameter('Steps', 50);
            p.addParameter('TVB', false);
            p.addParameter('TVD', false);
            p.addParameter('CFLMAX', 1.5);
            p.addParameter('CFLRefinement',0.1);
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
            
            %TODO: is this correct?
            assert(any([obj.testTVD obj.testTVB]),'Must choice TVD or TVB');
            
            
            %TODO : a better way to test for the problem
            obj.dudt = p.Results.integrator;
            obj.DT = obj.dudt.dfdx.dx*obj.cfl;
            [~, y] = obj.dudt.getState();
            obj.initTV = obj.calcTV(y);
            
            if obj.dudt.isSSP
                obj.sspCoef = obj.dudt.r;
            end
            
        end
        
        
        function run(obj, varargin)
            
            continueTest = 0;
            numViolation = 0;
            tvbSolution = [];
            tvdSolution = [];
            cfl_ = obj.cfl;
            dx_ = obj.dudt.dfdx.dx;
            while cfl_ < obj.CFLMAX
                obj.dudt.resetInitCondition();
                [~, y] = obj.dudt.getState();
                tv = obj.calcTV(y);
                tvmax = tv;
                dt = cfl_*dx_;
                ti = 0;
                tv_violation = false;
                TVMAX = [];
                
                for ti = 1:obj.Steps
                    obj.dudt.takeStep(dt);
                    [~, y] = obj.dudt.getState();
                    tv_new = obj.calcTV(y);
                    TVMAX = [TVMAX; tv_new];
                    tvmax = max(tv_new,tvmax);
                end
                
                if obj.testTVB
                    tvbSolution = [tvbSolution; [cfl_ max(TVMAX)]];
                end
                
                if obj.testTVD
                    tvdSolution = [tvdSolution; [cfl_ max(diff(TVMAX))]];
                end
                
                cfl_ = cfl_ + obj.cflRefine;
            end
            obj.TVB = tvbSolution;
            obj.TVD = tvdSolution;
        end
        
        function [ output ] = run_test(varargin) end
        
        function plotSolution(obj)
            
            %first calculate the SSP
            obj.calculateSSP();
            
            %now plot the solution
            semilogy(obj.TVD(:,1), obj.TVD(:,2),'x');
            ylim([1e-20 1]);
        end

        function calculateSSP(obj)
            if obj.testTVB
                diff_tv = log10(abs(obj.TVB(:,2) - obj.initTV));
                indTV = diff_tv < -14;
                diff_tv(indTV) = -15;
%                 badTV = diff_tv >= -5;
%                 diff_tv(badTV) = -5;
                indSSP = find(~indTV,1)-1;
                %plot(obj.TVB(:,1), diff_tv, 's', 'linewidth',2);
                obj.ssp = obj.TVB(indSSP,1);
            else
                goodIdx = obj.TVD(:,2) <= 1e-14;
                obj.TVD(goodIdx,2) = 1e-15;
                badIDX = obj.TVD(:,2) >= 1e-0;
                obj.TVD(badIDX,2) = 1e-0;
                indSSP = find(~goodIdx,1)-1;
                obj.ssp = obj.TVD(indSSP,1);
            end
        end
    end
    
    methods %(Access = protected)
        function log(obj, varargin) end
        
        function tv = calcTV(obj, u)
            tv = sum([abs([diff(u); (u(end)-u(1))])]);
        end
    end
    
end
