clear all; close all; clc


A = [0 0;1 0];
b = [1 0];
s = 2;
f = @(t,y) y;
y0 = 1;

dt = 0.1;
Tfinal = 1;
t = 0;

dudt = SSPTools.Integrators.ERK('A', A, 'b',b, 's', s, 'dydt', f, 'y0', y0);
dudx = SSPTools.Discretizers.FiniteDifference('bc','periodic');
dudx2 = SSPTools.Discretizers.Spectral('bc','periodic');

while t < Tfinal
    
    ynew = dudt.takeStep(dt);
    t = t+ dt;

    [t ynew]
end