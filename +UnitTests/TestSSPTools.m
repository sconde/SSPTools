classdef TestSSPTools < matlab.unittest.TestCase
    
    properties
        
    end
    
    % Test Method Block
    methods (Test)
        % includes unit test functions
        
        function testWENO5LinearAdvection(testCase)
            % test for completition
            
            import matlab.unittest.TestCase
            
            y0 = @(x) heaviside(x - (ceil((x+1)/2) -1)*2);
            imp_pro = TestProblems.PDEs.LinearAdvection('a',1);
            
            N = 100;
            dfdx = WenoCore.Weno5('N', N, 'domain', [-1, 1],...
                'kernel', 'WENO5', 'epsilon', 1e-16, 'p', 2,'Problem', imp_pro);
            
            dudt = SSPTools.Steppers.LoadERK('MethodName','FE',...
                'dfdx', dfdx,'y0', y0);
            
            dt = 0.1; t = 0;
            
            while t < 10*dt
                ynew = dudt.takeStep(dt);
                t = t+ dt;
            end
            
            testCase.assertTrue(true)
            assertTrue(testCase, true)
            
        end
        
        function testWENO5Burgers(testCase)
            % test for completition
            
            import matlab.unittest.TestCase
            
            y0 = @(x) heaviside(x - (ceil((x+1)/2) -1)*2);
            
            imp_pro = TestProblems.PDEs.Burgers();
            
            N = 100;
            dfdx = WenoCore.Weno5('N', N, 'domain', [-1, 1],...
                'kernel', 'WENO5', 'epsilon', 1e-16, 'p', 2,'Problem', imp_pro);
            
            dudt = SSPTools.Steppers.LoadERK('MethodName','FE',...
                'dfdx', dfdx,'y0', y0);
            
            dt = 0.1; t = 0;
            
            while t < 10*dt
                ynew = dudt.takeStep(dt);
                t = t+ dt;
            end
            
            testCase.assertTrue(true)
            assertTrue(testCase, true)
            
        end
        
        function testAdvectionDiffusion(testCase)
            
            
            y0 = @(x) sin(x);
            N = 8;
            exp_pro = TestProblems.PDEs.LinearAdvection('a',1);
            imp_pro = TestProblems.PDEs.LinearDiffusion();
            
            dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, 'N', N,...
                'Problem', exp_pro);
            
            
            dgdx = SSPTools.Discretizers.Spectral('derivativeOrder',2, 'N', N,...
                'Problem', imp_pro);
            
            sv = SSPTools.Steppers.LoadIMEX('MethodName','Stormer-Verlet',...
                'dfdx', dfdx,...
                'dgdx', dgdx, 'y0',y0);
            
            dt = 0.1; t = 0;
            
            while t < 10*dt
                ynew = sv.takeStep(dt);
                t = t+ dt;
            end
            
            testCase.assertTrue(true)
            assertTrue(testCase, true)
            
        end

    end
end
