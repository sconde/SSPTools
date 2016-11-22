%Sidafa Conde
%example of Mathieu ODE
clear all; close all; clc

addpath('../');

Tfinal = 10;

y0 = 1;

odePrb = TestProblems.ODEs.CurtissHirschfelder();

dudt = SSPTools.Steppers.LoadEmbeddedERK('MethodName','38rule43',...
    'ODE', odePrb, 'y0', y0, 'RelativeTolerance', 1e-4, 'AbsoluteTolerance', 1e-4,...
    'InitialStepSize', [],'VariableStepSize', true, 'Tfinal', Tfinal,...
    'FacMax',5);

T = []; DT = []; Y = []; ERR = []; badDT = [];
[t, y, dt, err, bad_dt] = dudt.getState();
T = [T; t]; DT = [DT; dt]; Y = [Y y]; ERR = [ERR; err]; badDT = [badDT; bad_dt];

while t < Tfinal
    
    dudt.takeStep(dt);
    [t, y, nextDt, err, bad_dt] = dudt.getState();
    dt = min(nextDt, Tfinal - t);
    T = [T; t]; DT = [DT; dt]; Y = [Y y]; ERR = [ERR; err]; badDT = [badDT; bad_dt];

end

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
% plotting the solution
subplot(3,1,1);
plot(T, Y(:,1), '-pk')
legend('y','location','southwest');
title('Solution')
xlim([0 Tfinal]);

% plotting the step size taken at each time step
subplot(3,1,2)
semilogy(T, DT,'-s')
hold on
semilogy(T, mean(DT)*ones(size(T)),'--k')
semilogy(badT, badDT, 'rx')
legend('dt','mean(dt)','rejected-dt','location','northwest');
title('Accepted Step Size')
ylim([0 5]);

% plotting the Error at each time-step
subplot(3,1,3)
semilogy(T, ERR,'-.k')
hold on
semilogy(T, 1e-4*ones(size(T)),'--')
legend('err','tolerance','location','northwest');
title('Error')
ylim([0 5]);
