classdef TestOrder < matlab.unittest.TestCase
    
    % Test Method Block
    methods (Test)
        % includes unit test functions
        
        function testBurgersRK2(testCase)
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            % test RK2 convegence for burgers
            dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, ...
                'N', 16);
            
            dudt = SSPTools.Steppers.LoadERK('MethodName','MidPoint',...
                'dfdx', dfdx, 'ExplicitProblem', TestProblems.PDEs.Burgers(),...
                'y0', @(x) sin(x));
            
            convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal',...
                0.4,'CFL', (1/2).^(1:5));
            
            convergencePDE.run();
            l2Error = convergencePDE.L2Error;
            DT = convergencePDE.DT(:);
            obsOrder = polyfit(log(DT),log(l2Error),1);
            obsOrder = obsOrder(1);
            
            testCase.assertThat(obsOrder, IsEqualTo(2, ...
                'Within', AbsoluteTolerance(0.2)))
        end
        
        function testBuckleeRK2(testCase)
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            % test RK2 convegence for BuckleyLeverett
            dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, ...
                'N', 16);
            
            dudt = SSPTools.Steppers.LoadERK('MethodName','MidPoint',...
                'dfdx', dfdx, 'ExplicitProblem', TestProblems.PDEs.BuckleyLeverett(),...
                'y0', @(x) sin(x));
            
            convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal',...
                0.4,'CFL', (1/2).^(1:5));
            
            convergencePDE.run();
            l2Error = convergencePDE.L2Error;
            DT = convergencePDE.DT(:);
            obsOrder = polyfit(log(DT),log(l2Error),1);
            obsOrder = obsOrder(1);
            
            testCase.assertThat(obsOrder, IsEqualTo(2, ...
                'Within', AbsoluteTolerance(0.2)))
        end
    end
end
