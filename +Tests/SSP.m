classdef SSP < Tests.Test
   
    properties
        cfl;
        t;
        dudt; % should be the integrator
        verbose;
        Tfinal;
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
            p.parse(varargin{:});
            
            %obj = obj@Tests.Test(varargin{:});
            
            obj.cfl = p.Results.cfl;
            obj.t = p.Results.t;
            obj.verbose = p.Results.verbose;
            obj.Tfinal = p.Results.Tfinal;
            obj.numViolation = p.Results.NumberAllowedViolation;
            
            %TODO : a better way to test for the problem
            obj.dudt = p.Results.integrator;
            obj.DT = obj.dudt.dfdx.dx*obj.cfl;
            [~, y] = obj.dudt.getState();
            obj.initTV = obj.calcTV(y);
        end
        
        
        function run(obj, varargin)
            
            continueTest = 0;
            numViolation = 0;
            tvSolution = [];
            cfl = obj.cfl;
            dx = obj.dudt.dfdx.dx;
            while numViolation < obj.numViolation
                obj.dudt.resetInitCondition();
                [~, y] = obj.dudt.getState();
                tv = obj.calcTV(y);
                dt = cfl*dx;
                ti = 0;
                tv_violation = false;
                TVMAX = [];
                
                for ti = 1:10
                    obj.dudt.takeStep(dt);
                    [~, y] = obj.dudt.getState();
                    tv_new = obj.calcTV(y);
                    TVMAX = [TVMAX; tv_new];
                end
                
                tvMax = max(TVMAX);
                tvSolution = [tvSolution; [cfl tvMax]];
                    
                if round(log10(tvMax- obj.initTV)) > -14
                    cfl = cfl - 0.001;
                    numViolation = numViolation + 1;
                else
                    cfl = cfl + 0.1;
                end
                
                continueTest = continueTest + 1;
            end
            obj.TV = tvSolution;
        end
        
        function [ output ] = run_test(varargin) end
        
        function ssp = plotSolution(obj)
            diff_tv = log10(abs(obj.TV(:,2) - obj.initTV));
            indTV = diff_tv < -14;
            diff_tv(indTV) = -15;
            indSSP = find(~indTV,1)-1;
            plot(obj.TV(:,1), diff_tv, 's', 'linewidth',2);
            ssp = obj.TV(indSSP,1);
        end
    end
    
    methods %(Access = protected)
        function log(obj, varargin) end 

        function tv = calcTV(obj, u)
              tv = sum([abs([diff(u); (u(end)-u(1))])]);
        end
    end
    
end
