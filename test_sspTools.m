clear all; close all; clc


A = [0 0;1 0];
b = [1 0];
s = 2;
N = 100;

dt = 1/N;
Tfinal = 1;
t = 0;

y0 = @(x) heaviside(x - (ceil((x+1)/2) -1)*2);
f = @(t, u) u;

prb = TestProblems.PDEs.PDE('f',f, 'y0', y0);
dudx = SSPTools.Discretizers.FiniteDifference('N', N, 'domain', [-1, 1],'bc','periodic');
dudt = SSPTools.Integrators.ERK('A', A, 'b',b, 's', s, 'dydt', f,...
    'dudx', dudx, 'problem', prb);

line1 = plot(dudx.x, dudt.y0,'-r','linewidth',2);
axis([-1 1 0 1]);

while t < Tfinal
    
    ynew = dudt.takeStep(dt);
    t = t+ dt;
    
    set(line1, 'ydata', ynew);
    drawnow;
    pause(0.1);
    
end