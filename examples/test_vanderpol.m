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
testing = 'DIRK';

y0 = [2; -0.6654321];

vdp = TestProblems.ODEs.Vanderpol();

if strcmpi(testing, 'erk')
    dudt = SSPTools.Steppers.ERK('A', A, 'b',b, 's', s,...
        'ExplicitProblem', vdp, 'y0', y0);
else
    dudt = SSPTools.Steppers.DIRK('A', At, 'b',bt, 's', s,...
        'ExplicitProblem', vdp, 'y0', y0);
end

plot(t, y0(1),'sr');
hold on
plot(t, y0(2), '*k');

while t < Tfinal
    
    y0 = dudt.takeStep(dt);
    t = t+ dt;
        
    plot(t, y0(1), 'sr');
    hold on
    plot(t, y0(2), '*k');
    pause(0.1);
end