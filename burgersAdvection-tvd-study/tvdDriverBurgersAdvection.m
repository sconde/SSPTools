clear all; close all; clc;

global dudt 

rk = loadMethod(2,2,2,2,1);
r = rk.r;
K = rk.k;

a = 1;
N = 601;

y0 = @(x) heaviside(x - (ceil((x+1)/2) -1)*2);


addpath('~/Dropbox/SSPTools/')

domain = [-1 1];

% define the advection flux, a = wave speed
imp_pro = TestProblems.PDEs.LinearAdvection('a',a);

% define the burgers flux
exp_pro = TestProblems.PDEs.Burgers();

% define the semi-discrete explcicit problem
dfdx = SSPTools.Discretizers.FiniteDifference('derivativeOrder',1, 'N', N,...
    'Problem', exp_pro, 'domain', domain);

% define the semi-discrete implicit problem
dgdx = SSPTools.Discretizers.FiniteDifference('derivativeOrder',1, 'N', N,...
    'Problem', imp_pro, 'domain', domain);

dudt = SSPTools.Steppers.IMEXRK('A', rk.A, 'b',rk.b, 'At', rk.At, 'bt', rk.bt,...
    'dfdx', dfdx,'dgdx', dgdx, 'y0', y0);


startingr = min(r,.005);
sp = min(r/10,.005);

lambda = sort(startingr + rand(1,20)*(r + 0.5 - startingr));
lambda = linspace(startingr, r + 0.5, 10);

cfl_refinement = max(diff(lambda))

tvdPDE = Tests.SSP('integrator', dudt,'r', r);

tvdPDE.run();
tvdPDE.plotSolution();
tvdPDE.ssp
