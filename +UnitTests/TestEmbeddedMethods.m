classdef TestEmbeddedMethods < matlab.unittest.TestCase
    
    properties
        
    end
    
    % Test Method Block
    methods (Test)
        % includes unit test functions
        
        function DormandPrince54(testCase)
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            t = 0;
            y0(1) = 0.994;
            y0(2) = 0.0;
            y0(3) = 0.0;
            y0(4) = -2.00158510637908252240537862224;
            Tfinal = 17.0652165601579625588917206249;
            
            aren = TestProblems.ODEs.Aren();
            
            dudt = SSPTools.Steppers.LoadEmbeddedERK('MethodName','DormandPrince54',...
                'ODE', aren, 'y0', y0, 'RelativeTolerance', 1e-7, 'AbsoluteTolerance', 1e-7,...
                'InitialStepSize', [],'VariableStepSize', true, 'Tfinal', Tfinal,'UseNew', false);
            [t, ~, dt, ~, ~] = dudt.getState();
            
            while t < Tfinal
                
                dudt.takeStep(dt);
                [t, ~, nextDt, ~, ~] = dudt.getState();
                
                dt = min(nextDt, Tfinal - t);
            end
            
              testCase.assertThat(dudt.nrejct, IsEqualTo(23)); % 22
              testCase.assertThat(dudt.naccpt, IsEqualTo(247)); %216
              testCase.assertThat(dudt.nfcn, IsEqualTo(1634)); %1442
        end
        
%         function testZonneveld43(testCase)
%             import matlab.unittest.TestCase
%             import matlab.unittest.constraints.IsEqualTo
%             import matlab.unittest.constraints.AbsoluteTolerance
%             
%             N = 16;
%             Tfinal = 0.4;
%             
%             
%             y0 = @(x) sin(x);
%             
%             
%             exp_pro = TestProblems.PDEs.LinearAdvection('a', 1);
%             dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, 'N', N,...
%                 'Problem', exp_pro);
%             
%             dudt = SSPTools.Steppers.LoadERK('MethodName', 'Zonneveld43',...
%                 'dfdx', dfdx, 'y0', y0, 'VariableStepSize', true);
%             
%             convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal', Tfinal,...
%                 'CFL', (1/2).^(1:5));
%             
%             espectedOrder = 4;
%             convergencePDE.run();
%             obsOrder = convergencePDE.getOrder('L2');
%             
%             testCase.verifyTrue(dudt.isEmbedded, true); % should be an embedded-RK
%             
%             testCase.assertThat(obsOrder, IsEqualTo(espectedOrder, ...
%                 'Within', AbsoluteTolerance(0.2)))
%         end
        
    end
end
