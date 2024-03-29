%Sidafa Conde
%example of vanderpol ODE
clear all; close all; clc


% TODO: fix this

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

%testing = 'ERK';
testing = 'DIRK';

y0 = 1;

vdp = TestProblems.ODEs.Dalquist('a', 2);

if strcmpi(testing, 'erk')
    dudt = SSPTools.Steppers.ERK('A', A, 'b',b, 's', s,...
        'ODE', vdp, 'y0', y0);
else
    dudt = SSPTools.Steppers.DIRK('A', At, 'b',bt, 's', s,...
        'ODE', vdp, 'y0', y0);
end

plot(t, y0(1),'sr');
hold on


while t < Tfinal
    
    y0 = dudt.takeStep(dt);
    t = t+ dt;
        
    plot(t, y0(1), 'sr');
    pause(0.1);
end
