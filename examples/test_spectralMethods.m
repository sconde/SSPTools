clear all; close all; clc

addpath('../')
A = [0 0;1 0]; b = [1 0]; s = 2; %TODO: infert s from size(A,1)
At = [0 0;0 1]; bt = [0 1];
N = 64;

dt = 0.01;
Tfinal = 2;
t = 0;

testing = 'ERK';

y0 = @(x) 0.9*sin(x); %heaviside(x - (ceil((x+1)/2) -1)*2);

imp_pro = TestProblems.PDEs.LinearAdvection('a', 1);

dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, 'N', N,...
    'Problem', imp_pro);

dudt = SSPTools.Steppers.ERK('A', A, 'b',b, 's', s,...
    'dfdx', dfdx, 'y0', y0);

line1 = plot(dfdx.x, dudt.y0(dfdx.x),'-r','linewidth',2);
axis([0 2*pi -1 1]);


while t < Tfinal
    
    ynew = dudt.takeStep(dt);
    t = t+ dt;
    
    set(line1, 'ydata', ynew);
    drawnow;
    pause(0.1);
    
end