clear all; close all; clc


A = [0 0;1 0]; b = [1 0]; s = 2; %TODO: infert s from size(A,1)
At = [0 0;0 1]; bt = [0 1];
N = 100;

dt = 0.001;
Tfinal = 2;
t = 0;

y0 = @(x) heaviside(x - (ceil((x+1)/2) -1)*2);
f = @(t, u) u;

imp_pro = TestProblems.PDEs.LinearAdvection('y0', y0, 'a', 1);
exp_pro = TestProblems.PDEs.Burgers('y0', y0);
exp_pro = TestProblems.PDEs.BuckleyLeverett('y0', y0);
dfdx = SSPTools.Discretizers.FiniteDifference('N', N, 'domain', [-1, 1],'bc','periodic');

dudt = SSPTools.Steppers.ERK('A', A, 'b',b, 's', s,...
    'dfdx', dfdx, 'ExplicitProblem', exp_pro);
% 
% dudt = SSPTools.Steppers.DIRK('A', A, 'b',b, 's', s,...
%     'dfdx', dfdx, 'ExplicitProblem', imp_pro);
% 
dudt = SSPTools.Steppers.IMEXRK('A', A, 'b',b, 's', s, 'At', At, 'bt', bt, 'dydt', f,...
    'dfdx', dfdx, 'ExplicitProblem', exp_pro, 'ImplicitProblem', imp_pro,...
    'dgdx', dfdx);

line1 = plot(dfdx.x, dudt.y0,'-r','linewidth',2);
axis([-1 1 0 1]);


while t < Tfinal
    
    ynew = dudt.takeStep(dt);
    t = t+ dt;
    
    set(line1, 'ydata', ynew);
    drawnow;
    pause(0.1);
    
end