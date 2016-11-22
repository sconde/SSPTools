clear all; close all; clc

% TODO: need to fix this

addpath('../')
A = [0 0;1 0]; b = [1 0]; s = 2; %TODO: infert s from size(A,1)
At = [0 0;0 1]; bt = [0 1];
N = 100;

dt = 0.001;
Tfinal = 2;
t = 0;

testing = 'euler';

if strcmpi(testing, 'euler')
    
    euler = TestProblems.PDEs.Euler1D('ProblemType', 'Sod','N',N);
    dfdx = WenoCore.Weno5('x',euler.x,...
        'kernel', 'WENO5', 'epsilon', 1e-16, 'p', 2,'Problem', euler);
    dudt = SSPTools.Steppers.LoadERK('MethodName','SSP54',...
    'dfdx', dfdx,'y0', euler.y0);
else
    
    euler = TestProblems.PDEs.Burgers();
    y0 = @(x) heaviside(x - (ceil((x+1)/2) -1)*2);
    
    dfdx = WenoCore.Weno5('N', N, 'domain', [-1, 1],...
        'kernel', 'WENO5', 'epsilon', 1e-16, 'p', 2,'Problem', euler);
    dudt = SSPTools.Steppers.LoadERK('MethodName','MidPoint',...
    'dfdx', dfdx,'y0',y0);
end

x = euler.x;
[~, Q] = dudt.getState();

t = 0;
subplot(3,1,1); rho_line = plot(x, Q(:,1), '-r'); %axis([0 1 0 1]);
subplot(3,1,2); vel_line = plot(x, Q(:,2), '-k'); %axis([0 1 0 0.4]);
subplot(3,1,3); p_line = plot(x, Q(:,3), '-b'); %axis([0 1 0 3]);
pause(0.1);

while t < Tfinal
    
    ynew = dudt.takeStep(dt);
    [~, Q] = dudt.getState();
    t = t+ dt;
    
    t
    set(rho_line, 'ydata', Q(:,1));
    set(vel_line, 'ydata', Q(:,2));
    set(p_line, 'ydata', Q(:,3));
    drawnow
    pause(0.1)
    pause(0.1);
    
end