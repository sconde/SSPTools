%Sidafa Conde
%example of vanderpol ODE
clear all; close all; clc

% TODO: fix ODEs

addpath('../');

Tfinal = 0.2;
N = 100;

euler = TestProblems.PDEs.Euler1D('ProblemType', 'Sod','N',N);

dfdx = WenoCore.Weno5('x',euler.x,...
    'kernel', 'WENO5', 'epsilon', 1e-16, 'p', 2,'Problem', euler);


dudt = SSPTools.Steppers.LoadEmbeddedERK('MethodName','38rule43',...
    'dfdx', dfdx,'y0', euler.y0,...
    'RelativeTolerance', 1e-4, 'AbsoluteTolerance', 1e-4,...
    'InitialStepSize', [],'VariableStepSize', true, 'Tfinal', Tfinal,...
    'FacMax',5);

T = []; DT = []; Y = []; ERR = []; badDT = [];
[t, y, dt, err, bad_dt] = dudt.getState();
T = [T; t]; DT = [DT; dt]; Y = [Y y]; ERR = [ERR; err]; badDT = [badDT; bad_dt];

x = euler.x;
subplot(3,1,1); rho_line = plot(x, y(:,1), '-r'); %axis([0 1 0 1]);
subplot(3,1,2); vel_line = plot(x, y(:,2), '-k'); %axis([0 1 0 0.4]);
subplot(3,1,3); p_line = plot(x, y(:,3), '-b'); %axis([0 1 0 3]);


while t < Tfinal
    
    dudt.takeStep(dt);
    [t, y, nextDt, err, bad_dt] = dudt.getState();
    dt = min(nextDt, Tfinal - t);
    T = [T; t]; DT = [DT; dt]; Y = [Y y]; ERR = [ERR; err]; badDT = [badDT; bad_dt];
    
    set(rho_line, 'ydata', y(:,1));
    set(vel_line, 'ydata', y(:,2));
    set(p_line, 'ydata', y(:,3));
    drawnow
    pause(0.1)

end

keyboard
close('all')
% print the summary
dudt.summary();

T = T(2:end); DT = DT(2:end); ERR = ERR(2:end); badDT = badDT(2:end);
% remove all the rejected Step and Solution

goodSol = isnan(badDT);
badDT = badDT(~goodSol);
BADERR = ERR(~goodSol);

[ismem, ind] = ismember(badDT,DT);
DT(ind) = [];
ERR(ind) = [];
Y = Y';
Y = Y(goodSol,:);
badT = T(~goodSol);
T = T(goodSol);

figure('Position',[100 100 1200 400]);

% plotting the step size taken at each time step
subplot(2,1,1)
semilogy(T, DT,'-s')
hold on
semilogy(T, mean(DT)*ones(size(T)),'--k')
semilogy(badT, badDT, 'rx')
legend('dt','mean(dt)','rejected-dt','location','northwest');
title('Accepted Step Size')
ylim([0 5]);

% plotting the Error at each time-step
subplot(2,1,2)
semilogy(T, ERR,'-.k')
hold on
semilogy(T, 1e-4*ones(size(T)),'--')
legend('err','tolerance','location','northwest');
title('Error')
ylim([0 5]);
