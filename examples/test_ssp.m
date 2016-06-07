clear all; close all; clc

addpath('../');
N = 1000;

y0 = @(x) heaviside(x - (ceil((x+1)/2) -1)*2);

exp_pro = TestProblems.PDEs.LinearAdvection('a', 1);

imp_pro = TestProblems.PDEs.LinearAdvection('a', 1);

dfdx = SSPTools.Discretizers.FiniteDifference('N', N, 'domain', [-1, 1],...
    'bc','periodic','Problem', exp_pro);

dgdx = SSPTools.Discretizers.FiniteDifference('N', N, 'domain', [-1, 1],...
    'bc','periodic','Problem', imp_pro);


dudt = SSPTools.Steppers.IMEXRK('A',0,'b',1,'At',1,'bt',1,...
    'dfdx', dfdx, 'dgdx', dgdx, 'y0', y0);

dudt.butcherCoef()

tvdPDE = Tests.SSP('integrator', dudt,'TVD',true,'CFLRefinement',0.1,...
    'CFLMAX',2,'CFL',1.0,'Steps',5);

tvdPDE.run();
tvdPDE.plotSolution();
tvdPDE.ssp
