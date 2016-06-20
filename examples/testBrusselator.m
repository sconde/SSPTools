%Sidafa Conde
%example of Brusselator ODE
% Purpose: recreate result from fig.4.1 in
% Solving Ordinary Differential Equation I
% Hairer. Pg 170

clear all; close all; clc
addpath('../');
addpath('~/Dropbox/SSPTools/');

dt = 0.01;
Tfinal = 20;
t = 0;

%method = 'DormandPrince54';
%method = 'Merson45';
%method = 'Zonneveld43';
%method = 'Felhberg45';
%method = '38rule43';
method = 'SSPEmbeddedRK';

%rk_method = '~/Dropbox/embedded-rk/butcher-optimization/Method/ERK/P4/Plin4/S5/method_typeG_r_0.6572570916924_acc_-16.mat';
rk_method = '~/Dropbox/embedded-rk/butcher-optimization/Method/ERK/P4/Plin5/S5/method_typeG_r_0.5271866339745_acc_-16.mat';
rk = load(rk_method);

y0 = [1.5; 3];

vdp = TestProblems.ODEs.Brusselator();


if strcmpi(method, 'sspembeddedrk')
    dudt = SSPTools.Steppers.EmbeddedERK('A', rk.A, 'b', rk.b, 'bhat', rk.bhat, 'p', rk.p, 'phat', rk.p-1,...
        'ODE', vdp, 'y0', y0, 'RelativeTolerance', 1e-4, 'AbsoluteTolerance', 1e-4,...
        'InitialStepSize', [],'VariableStepSize', true, 'Tfinal', Tfinal,...
        'FacMax',5);
else
    dudt = SSPTools.Steppers.LoadEmbeddedERK('MethodName',method,...
        'ODE', vdp, 'y0', y0, 'RelativeTolerance', 1e-4, 'AbsoluteTolerance', 1e-4,...
        'InitialStepSize', [],'VariableStepSize', true, 'Tfinal', Tfinal,...
        'FacMax',5);
end

T = []; DT = []; Y = []; ERR = []; badDT = [];

[t, y, dt, err, bad_dt] = dudt.getState();
T = [T; t]; DT = [DT; dt]; Y = [Y y]; ERR = [ERR; err]; badDT = [badDT; bad_dt];

while t < Tfinal
    
    dudt.takeStep(dt);
    [t, y, nextDt, err, bad_dt] = dudt.getState();
    dt = min(nextDt, Tfinal - t);
    T = [T; t]; DT = [DT; dt]; Y = [Y y]; ERR = [ERR; err]; badDT = [badDT; bad_dt];
end

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
title('Solution')
plot(T, Y(:,2), '-or')
hold on
plot(T, Y(:,1), '-pk')
xlim([0 Tfinal]);

% plotting the step size taken at each time step
subplot(3,1,2)
title('Accepted Step Size')
semilogy(T, DT,'-s')
hold on
semilogy(T, mean(DT)*ones(size(T)),'--k')
semilogy(badT, badDT, 'rx')
ylim([0 5]);

% plotting the Error at each time-step
subplot(3,1,3)
title('Error')
semilogy(T, ERR,'-.k')
hold on
semilogy(T, 1e-4*ones(size(T)),'--')
%semilogy(badT, BADERR, 'sr');
ylim([0 5]);
