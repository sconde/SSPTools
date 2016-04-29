clear all; close all; clc

addpath('../');
method = '~/Dropbox/imex-linear/src/butcher-optimization/Method/DIRK/G/Pex2/Pim2/Plin4/S5/K1/method_r1_0.5451888988915_acc_-15.mat';
rk = load(method);
A = rk.A; b = rk.b; s = rk.s;
At = rk.At; bt = rk.bt;

N = 100;

dt = 0.01;
Tfinal = 2;
t = 0;
f = @(u) u;

testing = 'ERK';
%testing = 'DIRK';
testing = 'IMEXRK';

%y0 = @(x) heaviside(x - (ceil((x+1)/2) -1)*2);
y0 = @(x) sin(x);

exp_pro = TestProblems.PDEs.LinearAdvection('a', 1);
imp_pro = TestProblems.PDEs.LinearDiffusion('nu', 0.1);

dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, 'N', N);
dgdx = SSPTools.Discretizers.Spectral('derivativeOrder',2, 'N', N);

dudt = SSPTools.Steppers.IMEXRK('A', A, 'b',b, 's', s, 'At', At, 'bt', bt, 'dydt', f,...
    'dfdx', dfdx, 'ExplicitProblem', exp_pro, 'ImplicitProblem', imp_pro,...
    'dgdx', dgdx, 'y0',y0);

line1 = plot(dfdx.x, dudt.y0(dfdx.x),'-r','linewidth',2);
axis([0 2*pi -1 1]);


while t < Tfinal
    
    ynew = dudt.takeStep(dt);
    t = t+ dt;
    
    set(line1, 'ydata', ynew);
    drawnow;
    pause(0.1);
    
end