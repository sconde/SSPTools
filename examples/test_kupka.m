%Sidafa Conde
% determine Error-Constant
% Total-Variation-Diminishing Implicit-Explicit Runge-Kutta Methods for the
% Simulation of Double-Diffusive Convection in Astrophyics (Friedrich
% Kupka)

clear all; close all; clc

% TODO: fix this

addpath('../');
method = '~/Dropbox/imex-linear/src/butcher-optimization/Method/DIRK/G/Pex2/Pim2/Plin4/S5/K1/method_r1_0.5451888988915_acc_-15.mat';
rk = load(method);
A = rk.A; b = rk.b; s = rk.s;
At = rk.At; bt = rk.bt;

N = 100;

dt = 0.01;
Tfinal = 1.3;
t = 0;

y0 = 0;

exp_pro = TestProblems.ODEs.KupkaExplicit();
imp_pro = TestProblems.ODEs.KupkaImplicit();


dudt = SSPTools.Steppers.LoadIMEX('MethodName', 'IMEXSSP2222',...
    'ODE', exp_pro, 'ImplicitODE', imp_pro, 'y0', y0);


% convergence test
convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal', Tfinal,...
    'DT', Tfinal./[40:10:100],'ExactSolution',@(t) tan(t));
convergencePDE.run();
convergencePDE.complete();
l2Error = convergencePDE.getError('l2');
DT = convergencePDE.getDT();

ErrorConst = max(l2Error(1:end-1)./l2Error(2:end));
l2Error./DT(1)
% hfig = gcf;
% axesObjs = get(hfig, 'Children');  %axes handles
% dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
% objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
% xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
% xt = cell2mat(xdata);
% plot(xt, tan(xt), '-k');
