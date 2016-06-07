classdef TestODEs < matlab.unittest.TestCase
    
    properties
        
    end
    
    % Test Method Block
    methods (Test)
        % includes unit test functions
        
        function testDalquitsExplicit(testCase)
            % test for completition
            
            import matlab.unittest.TestCase
            
            dt = 0.01;
            Tfinal = 2;
            t = 0;
            
            y0 = 1;
            
            dalq = TestProblems.ODEs.Dalquist('a', 2);
            
            dudt = SSPTools.Steppers.LoadERK('MethodName','Heuns',...
                'ODE', dalq, 'y0', y0);
            
            while t < Tfinal
                
                dudt.takeStep(dt);
                t = t+ dt;
            end
            
            testCase.assertTrue(true)
            assertTrue(testCase, true)
            
        end
        
        function testVanderpolDIRK(testCase)
            % test for completition
            
            import matlab.unittest.TestCase
            
            dt = 0.01;
            Tfinal = 2;
            t = 0;
            
            y0 = [2; -0.6654321];
            
            vdp = TestProblems.ODEs.Vanderpol();
            
            dudt = SSPTools.Steppers.DIRK('A',[0 0;0 1], 'b', [1/2 1/2],...
                'ODE', vdp, 'y0', y0);
            
            while t < Tfinal
                
                dudt.takeStep(dt);
                t = t+ dt;
            end
            
            testCase.assertTrue(true);
            assertTrue(testCase, true);
            
        end
        
        function testVanderpolSDIRK23(testCase)
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            Tfinal = 2;

            y0 = [2; -0.6654321];
            
            vdp = TestProblems.ODEs.Vanderpol();
            
            dudt = SSPTools.Steppers.LoadDIRK('MethodName','SDIRK23',...
                'ODE', vdp, 'y0', y0);
            
            convergence_pro = Tests.Convergence('integrator', dudt,'Tfinal', Tfinal,...
                'DT', Tfinal./(200:20:300));
            
            convergence_pro.run();
            obsOrder = convergence_pro.getOrder('l2');
            expectedOrder = 3;
            testCase.assertThat(obsOrder, IsEqualTo(expectedOrder, ...
                'Within', AbsoluteTolerance(0.2)))
        end
        
        function testDalquitsFelhberg45(testCase)
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            Tfinal = 1;
            
            y0 = [2; -0.6654321];
            
            vdp = TestProblems.ODEs.Vanderpol();
            
            dudt = SSPTools.Steppers.LoadERK('MethodName','Felhberg45',...
                'ODE', vdp, 'y0', y0);
            
            convergence_pro = Tests.Convergence('integrator', dudt,'Tfinal', Tfinal,...
                'DT', Tfinal./(160:20:220));
            
            convergence_pro.run();
            obsOrder = convergence_pro.getOrder('l2');
            expectedOrder = 4;
            testCase.assertThat(obsOrder, IsEqualTo(expectedOrder, ...
                'Within', AbsoluteTolerance(0.2)))
        end
        
          function testVanderpolSDIRK22(testCase)
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            Tfinal = 2;
            
            y0 = [2; -0.6654321];
            
            vdp = TestProblems.ODEs.Vanderpol();
            
            dudt = SSPTools.Steppers.LoadDIRK('MethodName','SDIRK22',...
                'ODE', vdp, 'y0', y0);
            
            convergence_pro = Tests.Convergence('integrator', dudt,'Tfinal', Tfinal,...
                'DT', Tfinal./(160:20:220));
            
            convergence_pro.run();
            obsOrder = convergence_pro.getOrder('l2');
            expectedOrder = 2;
            testCase.assertThat(obsOrder, IsEqualTo(expectedOrder, ...
                'Within', AbsoluteTolerance(0.2)))
        end
        
        function testKupkaOrderSSP22(testCase)
            % ODE IMEX problem
            
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            y0 = 0;
            Tfinal = 1.3;
            
            exp_pro = TestProblems.ODEs.KupkaExplicit();
            imp_pro = TestProblems.ODEs.KupkaImplicit();
            
            dudt = SSPTools.Steppers.LoadIMEX('MethodName', 'IMEXSSP2222',...
                'ODE', exp_pro, 'ImplicitODE', imp_pro, 'y0', y0);
            
            convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal', Tfinal,...
                'DT', Tfinal./(40:10:100),'ExactSolution',@(t) tan(t));
            convergencePDE.run();
            obsOrder = convergencePDE.getOrder('l2');
            testCase.assertThat(obsOrder, IsEqualTo(2, ...
                'Within', AbsoluteTolerance(0.2)))
        end
        
    end
end
