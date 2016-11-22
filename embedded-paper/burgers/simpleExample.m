clear all; close all; clc

% TODO: need to fix this

addpath('../')
A = [0 0;1 0]; b = [1 0]; s = 2; %TODO: infert s from size(A,1)
At = [0 0;0 1]; bt = [0 1];
N = 100;

dt = 0.01;
Tfinal = 2;
t = 0;

testing = 'ERK';

y0 = @(x) heaviside(x - (ceil((x+1)/2) -1)*2);

imp_pro = TestProblems.PDEs.Burgers();

dfdx = WenoCore.Weno5('N', N, 'domain', [-1, 1],...
    'kernel', 'WENO5', 'epsilon', 1e-16, 'p', 2,'Problem', imp_pro);

dudt = SSPTools.Steppers.LoadERK('MethodName','MidPoint',...
    'dfdx', dfdx,'y0', y0);

line1 = plot(dfdx.x, dudt.y0(dfdx.x),'-r','linewidth',2);
%axis([-1 1 0 1]);

%keyboard

while t < Tfinal
    
    ynew = dudt.takeStep(dt);
    t = t+ dt;
    
    set(line1, 'ydata', ynew);
    drawnow;
    pause(0.1);
    
end