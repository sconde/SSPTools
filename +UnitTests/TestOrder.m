classdef TestOrder < matlab.unittest.TestCase
    
    properties
        
    end
    
    % Test Method Block
    methods (Test)
        % includes unit test functions
        
        function testBurgersRK2(testCase)
            % test RK2 convegence for burgers
            
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, ...
                'N', 16, 'Problem',TestProblems.PDEs.Burgers());
            
            dudt = SSPTools.Steppers.LoadERK('MethodName','MidPoint',...
                'dfdx', dfdx,...
                'y0', @(x) sin(x));
            
            convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal',...
                0.4,'CFL', (1/2).^(1:5));
            
            convergencePDE.run();
            obsOrder = convergencePDE.getOrder('L2');
            
            testCase.assertThat(obsOrder, IsEqualTo(2, ...
                'Within', AbsoluteTolerance(0.2)))
        end
        
        function testBuckleeRK2(testCase)
            % test RK2 convegence for BuckleyLeverett
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, ...
                'N', 16,'Problem',TestProblems.PDEs.BuckleyLeverett());
            
            dudt = SSPTools.Steppers.LoadERK('MethodName','MidPoint',...
                'dfdx', dfdx,...
                'y0', @(x) sin(x));
            
            convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal',...
                0.4,'CFL', (1/2).^(1:5));
            
            convergencePDE.run();
            obsOrder = convergencePDE.getOrder('L2');
            
            testCase.assertThat(obsOrder, IsEqualTo(2, ...
                'Within', AbsoluteTolerance(0.2)))
        end
        
        function testAdvectionDiffusion(testCase)
            % test advection diffusion with ssp33
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            y0 = @(x) sin(x);
            N = 8;
            Tfinal = 0.4;
            
            exp_pro = TestProblems.PDEs.LinearAdvection('a', 1);
            imp_pro = TestProblems.PDEs.LinearDiffusion('nu', 0.1);
            
            dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, 'N', N);
            
            dgdx = SSPTools.Discretizers.Spectral('derivativeOrder',2, 'N', N);
            
            
            CFL = (1/2).^(1:5);
            
            
            ssp22 = SSPTools.Steppers.LoadIMEX('MethodName','IMEXSSP3333',...
                'dfdx', dfdx,'ExplicitProblem', exp_pro,...
                'dgdx', dgdx, 'y0',y0,'ImplicitProblem', imp_pro);
            
            
            ssp22_convergence = Tests.Convergence('integrator', ssp22,'Tfinal', Tfinal,...
                'CFL', CFL);
            
            
            ssp22_convergence.run();
            obsOrder = ssp22_convergence.getOrder('L2');
            testCase.assertThat(obsOrder, IsEqualTo(3, ...
                'Within', AbsoluteTolerance(0.2)))
        end
    end
end
