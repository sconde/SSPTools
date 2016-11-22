%Sidafa Conde
%example of Brusselator ODE

clear all; close all; clc
addpath('../');

Tfinal = 10;

method = 'DormandPrince54';

y0 = [1; 0];

vdp = TestProblems.ODEs.Mathieu();

p = 4; phat = 5;
A = [0 0 0 0 0;
    1/3 0 0 0 0;
    1/6 1/6 0 0 0;
    1/8 0 3/8 0 0;
    1/2 0 -3/2 2 0];
b = [1/6 0  0 2/3 1/6];
bhat = [1/10 0 3/10 2/5 1/5];


dudt = SSPTools.Steppers.EmbeddedERK('A',A,'b',b,'bhat',bhat, 'p', p, 'phat', phat,...
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
