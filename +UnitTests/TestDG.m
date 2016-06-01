classdef TestDG < matlab.unittest.TestCase
    
    properties
        
    end
    
    % Test Method Block
    methods (Test)
        % includes unit test functions
        
        function testDGLinearAdvection(testCase)
            % test for completition
            
            import matlab.unittest.TestCase
            
            dt = 0.01;
            Tfinal = 10*dt;
            
            y0 = @(x) sin(pi*x);
            
            nE = 20;    % number of elements
            K  = 2;     % degree of accuracy
            
            exp_pro = TestProblems.PDEs.LinearAdvection('a',1);
            
            xmesh = NDG.Mesh1D('Domain', [-1 1], 'NumberElements', nE,...
                'SolutionDegree', K, 'QuadratureType', 'LGL');
            
            
            dg = NDG.NDG('Mesh', xmesh, 'Problem', exp_pro);
            
            dudt = SSPTools.Steppers.LoadERK('MethodName', 'FE',...
                'dfdx', dg, 'y0', y0);
            
            t = 0;
            while t < Tfinal
                ynew = dudt.takeStep(dt);
                t = t+ dt;
            end
            
            testCase.assertTrue(true)
            assertTrue(testCase, true)
            
        end
        
        function testDGBurgers(testCase)
            % test for completition
            
            import matlab.unittest.TestCase
            
            dt = 0.01;
            Tfinal = 10*dt;
            
            y0 = @(x) sin(pi*x);
            
            nE = 20;    % number of elements
            K  = 2;     % degree of accuracy
            
            exp_pro = TestProblems.PDEs.Burgers();
            
            xmesh = NDG.Mesh1D('Domain', [-1 1], 'NumberElements', nE,...
                'SolutionDegree', K, 'QuadratureType', 'LGL');
            
            
            dg = NDG.NDG('Mesh', xmesh, 'Problem', exp_pro);
            
            dudt = SSPTools.Steppers.LoadERK('MethodName', 'FE',...
                'dfdx', dg, 'y0', y0);
 
            t = 0;
            while t < Tfinal
                ynew = dudt.takeStep(dt);
                t = t+ dt;
            end
            
            testCase.assertTrue(true)
            assertTrue(testCase, true)
            
        end
        
       
    end
end
