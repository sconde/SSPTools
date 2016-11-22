%Sidafa Conde
%example of Brusselator ODE

clear all; close all; clc
addpath('../');

Tfinal = 1;

method = 'DormandPrince54';
N = 100;

y0 = @(x) heaviside(x - (ceil((x+1)/2) -1)*2);

imp_pro = TestProblems.PDEs.LinearAdvection('a',1);

dfdx = WenoCore.Weno5('N', N, 'domain', [-1, 1],...
    'kernel', 'WENO5', 'epsilon', 1e-16, 'p', 2,'Problem', imp_pro);

dudt = SSPTools.Steppers.LoadEmbeddedERK('MethodName',method,...
    'dfdx', dfdx,'y0', y0,...
    'RelativeTolerance', 1e-4, 'AbsoluteTolerance', 1e-4,...
    'InitialStepSize', [],'VariableStepSize', true, 'Tfinal', Tfinal,...
    'FacMax',5,'UseNew', false);

toleranceRange = logspace(-1,-8,5);
embeddedTest = Tests.EmbeddedRK('integrator', dudt,'Tolerance',toleranceRange);

embeddedTest.run();

numberAccp = embeddedTest.NAccept;
numberReject = embeddedTest.NReject;
numberFncEval = embeddedTest.NFunEval;
numberSteps = embeddedTest.NStep;
tolerance = embeddedTest.TOL;
info = [tolerance' numberFncEval' numberSteps' numberAccp' numberReject'];

fprintf(1,'%10s %10s %10s %10s %10s\n','tolerance', 'FncEval', 'Steps',...
    'Accepted','Rejected');
fprintf(1,'%10.2e %10d %10d %10d %10d\n',info');
