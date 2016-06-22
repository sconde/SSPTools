classdef TestOrder < matlab.unittest.TestCase
    
    properties
        
    end
    
    % Test Method Block
    methods (Test)
        % includes unit test functions
        function testBurgersSDIRK22(testCase)
            % test RK2 convegence for burgers
            
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, ...
                'N', 16, 'Problem',TestProblems.PDEs.Burgers());
            
            dudt = SSPTools.Steppers.LoadDIRK('MethodName','SDIRK22',...
                'dfdx', dfdx,...
                'y0', @(x) sin(x));
            
            convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal',...
                0.4,'CFL', (1/2).^(1:5));
            
            convergencePDE.run();
            obsOrder = convergencePDE.getOrder('L2');
            espectedOrder = 2;
            
            testCase.assertThat(obsOrder, IsEqualTo(espectedOrder, ...
                'Within', AbsoluteTolerance(0.2)))
        end
        
        function testBurgersSDIRK23(testCase)
            % test RK2 convegence for burgers
            
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, ...
                'N', 16, 'Problem',TestProblems.PDEs.Burgers());
            
            dudt = SSPTools.Steppers.LoadDIRK('MethodName','SDIRK23',...
                'dfdx', dfdx,...
                'y0', @(x) sin(x));
            
            convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal',...
                0.4,'CFL', (1/2).^(1:5));
            
            convergencePDE.run();
            obsOrder = convergencePDE.getOrder('L2');
            espectedOrder = 3;
            
            testCase.assertThat(obsOrder, IsEqualTo(espectedOrder, ...
                'Within', AbsoluteTolerance(0.2)))
        end
        
%         function testBurgersSDIRK54(testCase) % TODO: not correct
%             % test RK2 convegence for burgers
%             
%             import matlab.unittest.TestCase
%             import matlab.unittest.constraints.IsEqualTo
%             import matlab.unittest.constraints.AbsoluteTolerance
%             
%             dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, ...
%                 'N', 16, 'Problem',TestProblems.PDEs.Burgers());
%             
%             dudt = SSPTools.Steppers.LoadDIRK('MethodName','SDIRK54',...
%                 'dfdx', dfdx,...
%                 'y0', @(x) sin(x));
%             
%             convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal',...
%                 0.4,'CFL', (1/2).^(1:5));
%             
%             convergencePDE.run();
%             obsOrder = convergencePDE.getOrder('L2');
%             espectedOrder = 4;
%             
%             testCase.assertThat(obsOrder, IsEqualTo(espectedOrder, ...
%                 'Within', AbsoluteTolerance(0.2)))
%         end
        
        function testBurgersSDIRK34(testCase)
            % test RK2 convegence for burgers
            
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, ...
                'N', 16, 'Problem',TestProblems.PDEs.Burgers());
            
            dudt = SSPTools.Steppers.LoadDIRK('MethodName','SDIRK34',...
                'dfdx', dfdx,...
                'y0', @(x) sin(x));
            
            convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal',...
                0.4,'CFL', (1/2).^(1:5));
            
            convergencePDE.run();
            obsOrder = convergencePDE.getOrder('L2');
            espectedOrder = 4;
            
            testCase.assertThat(obsOrder, IsEqualTo(espectedOrder, ...
                'Within', AbsoluteTolerance(0.2)))
        end
        
        function testBurgersSSP54(testCase)
            % test RK2 convegence for burgers
            
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, ...
                'N', 16, 'Problem',TestProblems.PDEs.Burgers());
            
            dudt = SSPTools.Steppers.LoadERK('MethodName','SSP54',...
                'dfdx', dfdx,...
                'y0', @(x) sin(x));
            
            convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal',...
                0.4,'CFL', (1/2).^(1:5));
            
            convergencePDE.run();
            obsOrder = convergencePDE.getOrder('L2');
            espectedOrder = 4;
            
            testCase.assertThat(obsOrder, IsEqualTo(espectedOrder, ...
                'Within', AbsoluteTolerance(0.2)))
        end
        
        function testBurgersSSP33(testCase)
            % test RK2 convegence for burgers
            
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, ...
                'N', 16, 'Problem',TestProblems.PDEs.Burgers());
            
            dudt = SSPTools.Steppers.LoadERK('MethodName','SSP33',...
                'dfdx', dfdx,...
                'y0', @(x) sin(x));
            
            convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal',...
                0.4,'CFL', (1/2).^(1:5));
            
            convergencePDE.run();
            obsOrder = convergencePDE.getOrder('L2');
            espectedOrder = 3;
            
            testCase.assertThat(obsOrder, IsEqualTo(espectedOrder, ...
                'Within', AbsoluteTolerance(0.2)))
        end
        
        function testBurgersSSP22(testCase)
            % test RK2 convegence for burgers
            
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, ...
                'N', 16, 'Problem',TestProblems.PDEs.Burgers());
            
            dudt = SSPTools.Steppers.LoadERK('MethodName','SSP22',...
                'dfdx', dfdx,...
                'y0', @(x) sin(x));
            
            convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal',...
                0.4,'CFL', (1/2).^(1:5));
            
            convergencePDE.run();
            obsOrder = convergencePDE.getOrder('L2');
            espectedOrder = 2;
            
            testCase.assertThat(obsOrder, IsEqualTo(espectedOrder, ...
                'Within', AbsoluteTolerance(0.2)))
        end
        
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
            
            dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, 'N', N,...
                'Problem', exp_pro);
            
            dgdx = SSPTools.Discretizers.Spectral('derivativeOrder',2, 'N', N,...
                'Problem', imp_pro);
            
            
            CFL = (1/2).^(1:5);
            
            
            ssp22 = SSPTools.Steppers.LoadIMEX('MethodName','IMEXSSP3333',...
                'dfdx', dfdx,...
                'dgdx', dgdx, 'y0',y0);
            
            
            ssp22_convergence = Tests.Convergence('integrator', ssp22,'Tfinal', Tfinal,...
                'CFL', CFL);
            
            
            ssp22_convergence.run();
            obsOrder = ssp22_convergence.getOrder('L2');
            testCase.assertThat(obsOrder, IsEqualTo(3, ...
                'Within', AbsoluteTolerance(0.2)))
        end
        
        function testBurgersRK4(testCase)
            % test RK2 convegence for burgers
            
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, ...
                'N', 16, 'Problem',TestProblems.PDEs.Burgers());
            
            dudt = SSPTools.Steppers.LoadERK('MethodName','RK4',...
                'dfdx', dfdx,...
                'y0', @(x) sin(x));
            
            convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal',...
                0.4,'CFL', (1/2).^(1:5));
            
            convergencePDE.run();
            obsOrder = convergencePDE.getOrder('L2');
            
            testCase.assertThat(obsOrder, IsEqualTo(4, ...
                'Within', AbsoluteTolerance(0.2)))
        end
    end
end
