clear all; close all; clc

addpath('../');
N = 1000;

y0 = @(x) heaviside(x - (ceil((x+1)/2) -1)*2);

exp_pro = TestProblems.PDEs.LinearAdvection('a', 1);

imp_pro = TestProblems.PDEs.LinearAdvection('a', 10);

dfdx = SSPTools.Discretizers.FiniteDifference('N', N, 'domain', [-1, 1],...
    'bc','periodic','Problem', exp_pro);

dgdx = SSPTools.Discretizers.FiniteDifference('N', N, 'domain', [-1, 1],...
    'bc','periodic','Problem', imp_pro);


dudt = SSPTools.Steppers.LoadERK('MethodName', 'FE',...
    'dfdx', dfdx, 'dgdx', dgdx, 'y0', y0);

dudt.butcherCoef()

tvdPDE = Tests.SSP('integrator', dudt,'TVD',true,'CFLRefinement',0.01,...
    'CFLMAX',1.1,'CFL',0.01);

tvdPDE.run();
tvdPDE.plotSolution();
tvdPDE.ssp
