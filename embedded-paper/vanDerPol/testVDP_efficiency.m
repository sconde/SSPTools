%Sidafa Conde
%example of Brusselator ODE

clear all; close all; clc
addpath('../');

Tfinal = 2;

method = 'DormandPrince54';

y0 = [2; -0.6654321];

vdp = TestProblems.ODEs.Vanderpol();


dudt = SSPTools.Steppers.LoadEmbeddedERK('MethodName',method,...
    'ODE', vdp, 'y0', y0, 'RelativeTolerance', 1e-4, 'AbsoluteTolerance', 1e-4,...
    'InitialStepSize', [],'VariableStepSize', true, 'Tfinal', Tfinal,...
    'FacMax',5);

toleranceRange = logspace(-1,-8,10);
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
