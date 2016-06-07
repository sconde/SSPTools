classdef TestEmbeddedMethods < matlab.unittest.TestCase
    
    properties
        
    end
    
    % Test Method Block
    methods (Test)
        % includes unit test functions
        
        function testMerson45(testCase)
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            N = 16;
            Tfinal = 0.4;
            
            
            y0 = @(x) sin(x);
            
            
            exp_pro = TestProblems.PDEs.LinearAdvection('a', 1);
            dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, 'N', N,...
                'Problem', exp_pro);
            
            dudt = SSPTools.Steppers.LoadERK('MethodName', 'Merson45',...
                'dfdx', dfdx, 'y0', y0, 'VariableStepSize', true);
            
            convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal', Tfinal,...
                'CFL', (1/2).^(1:5));
            
            espectedOrder = 4;
            convergencePDE.run();
            obsOrder = convergencePDE.getOrder('L2');
            
            testCase.verifyTrue(dudt.isEmbedded, true); % should be an embedded-RK
            
            testCase.assertThat(obsOrder, IsEqualTo(espectedOrder, ...
                'Within', AbsoluteTolerance(0.2)))
        end
        
        function testZonneveld43(testCase)
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            N = 16;
            Tfinal = 0.4;
            
            
            y0 = @(x) sin(x);
            
            
            exp_pro = TestProblems.PDEs.LinearAdvection('a', 1);
            dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, 'N', N,...
                'Problem', exp_pro);
            
            dudt = SSPTools.Steppers.LoadERK('MethodName', 'Zonneveld43',...
                'dfdx', dfdx, 'y0', y0, 'VariableStepSize', true);
            
            convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal', Tfinal,...
                'CFL', (1/2).^(1:5));
            
            espectedOrder = 4;
            convergencePDE.run();
            obsOrder = convergencePDE.getOrder('L2');
            
            testCase.verifyTrue(dudt.isEmbedded, true); % should be an embedded-RK
            
            testCase.assertThat(obsOrder, IsEqualTo(espectedOrder, ...
                'Within', AbsoluteTolerance(0.2)))
        end
        
    end
end
