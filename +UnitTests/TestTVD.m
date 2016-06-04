classdef TestTVD < matlab.unittest.TestCase
    
    properties
        
    end
    
    % Test Method Block
    methods (Test)
        % includes unit test functions
        
        function testSSP22(testCase)
            % test SSP22 SSP for square wave
            
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            y0 = @(x) heaviside(x - (ceil((x+1)/2) -1)*2);

            exp_pro = TestProblems.PDEs.LinearAdvection('a', 1);
            
            dfdx = SSPTools.Discretizers.FiniteDifference('N', 500,...
                'domain', [-1, 1],'bc','periodic','Problem', exp_pro);
            
            dudt = SSPTools.Steppers.LoadERK('MethodName', 'SSP22',...
                'dfdx', dfdx,  'y0', y0);
            
            tvdPDE = Tests.SSP('integrator', dudt,'TVD',true,'CFLRefinement',0.001,...
                'CFLMAX',1.06,'CFL',0.85);
            
            tvdPDE.run();
            tvdPDE.calculateSSP();
            
            testCase.assertThat(tvdPDE.ssp, IsEqualTo(1, ...
                'Within', AbsoluteTolerance(0.2)))
        end
        
        function testFE(testCase)
            % test FE SSP for square wave
            
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            y0 = @(x) heaviside(x - (ceil((x+1)/2) -1)*2);

            exp_pro = TestProblems.PDEs.LinearAdvection('a', 1);
            
            dfdx = SSPTools.Discretizers.FiniteDifference('N', 500,...
                'domain', [-1, 1],'bc','periodic','Problem', exp_pro);
            
            dudt = SSPTools.Steppers.LoadERK('MethodName', 'FE',...
                'dfdx', dfdx,  'y0', y0);
            
            tvdPDE = Tests.SSP('integrator', dudt,'TVD',true,'CFLRefinement',0.001,...
                'CFLMAX',1.06,'CFL',0.85);
            
            tvdPDE.run();
            tvdPDE.calculateSSP();
            
            testCase.assertThat(tvdPDE.ssp, IsEqualTo(1, ...
                'Within', AbsoluteTolerance(0.2)))
        end
        
        
    end
end
