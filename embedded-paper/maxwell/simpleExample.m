clear all; close all; clc

addpath('../');


% Set Problem paramters
L = 1;
FinalTime = 0.2;
N = 100;

dt = 0.01;
Tfinal = 1;

maxwell = TestProblems.PDEs.Maxwell1D('ProblemType', 'Cavity','N',N);

dfdx = WenoCore.Weno5('x',maxwell.x,...
    'kernel', 'WENO5', 'epsilon', 1e-16, 'p', 2,'Problem', maxwell);

dudt = SSPTools.Steppers.LoadERK('MethodName','SSP54',...
    'dfdx', dfdx,'y0', maxwell.y0);
x = maxwell.x;
[~, Q] = dudt.getState();

t = 0;
subplot(2,1,1); rho_line = plot(x, Q(:,1), '-r'); %axis([0 1 0 1]);
subplot(2,1,2); vel_line = plot(x, Q(:,2), '-k'); %axis([0 1 0 0.4]);

pause(0.1);

while t < Tfinal
    dudt.takeStep(dt);
    t = t+ dt;
    [~, Q] = dudt.getState();
    set(rho_line, 'ydata', Q(:,1));
    set(vel_line, 'ydata', Q(:,2));
    drawnow
    pause(0.1)
end


