%Sidafa Conde
%example of vanderpol ODE
clear all; close all; clc


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


dudt = SSPTools.Steppers.LoadIMEX('MethodName', 'IMEXSSP1111LPM',...
    'ExplicitProblem', exp_pro, 'ImplicitProblem', imp_pro, 'y0', y0);

% plot(t, y0(1),'sr');
% hold on
% 
% while t < Tfinal
%     
%     y0 = dudt.takeStep(dt);
%     t = t+ dt;
%         
%     plot(t, y0(1), 'sr');
%     pause(0.1);
% end

% convergence test
convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal', Tfinal,...
    'DT', Tfinal./[40:10:80],'ExactSolution',@(t) tan(t));
convergencePDE.run();
convergencePDE.complete()
l2Error = convergencePDE.getError('l2')

% hfig = gcf;
% axesObjs = get(hfig, 'Children');  %axes handles
% dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
% objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
% xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
% xt = cell2mat(xdata);
% plot(xt, tan(xt), '-k');
