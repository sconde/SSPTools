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

% method = '38Rule43';                % 106 accepted + 20 rejected
% method = 'DormandPrince54';
% method = 'Felhberg45';                 % 62 accepted + 39 rejected
% method = 'Zonneveld43';
% method = 'Merson45';                   % 77 accepted + 20 rejected
% method = 'HeunEuler';                  % 1632 accepted + 0 rejected
% method = 'BogackiShampine';            % 166 accepted + 36 rejected
method = 'SSPEmbedded';

%rk = '~/Documents/embedded-rk/butcher-optimization/Method/ERK/P4/Plin4/PhatLn4/S5/method_typeG_r_1.5061659316780_acc_-16.mat';
%rk = load(rk);

y0 = [1.5; 3];

vdp = TestProblems.ODEs.Brusselator();

% if strcmpi(method, '38Rule43')
%     
%     dudt = SSPTools.Steppers.LoadEmbeddedERK('MethodName', '38Rule43',...
%         'ODE', vdp, 'y0', y0, 'RelTol', 1e-4, 'AbsTol', 1e-4,...
%         'InitialStepSize', [],'VariableStepSize', true,'MaxStepSize', 1e1);
%     
% elseif strcmpi(method, 'DormandPrince54')
%     dudt = SSPTools.Steppers.LoadEmbeddedERK('MethodName', 'DormandPrince54',...
%         'ODE', vdp, 'y0', y0, 'RelTol', 1e-4, 'AbsTol', 1e-4,...
%         'InitialStepSize', [],'VariableStepSize', true,'MaxStepSize', 1e1);
% elseif strcmpi(method, 'Felhberg45')
%     dudt = SSPTools.Steppers.LoadEmbeddedERK('MethodName', 'Felhberg45',...
%         'ODE', vdp, 'y0', y0, 'RelTol', 1e-4, 'AbsTol', 1e-4,...
%         'InitialStepSize', [],'VariableStepSize', true,'MaxStepSize', 1e1);
% elseif strcmpi(method, 'Zonneveld43')
%     dudt = SSPTools.Steppers.LoadEmbeddedERK('MethodName', 'Zonneveld43',...
%         'ODE', vdp, 'y0', y0, 'RelTol', 1e-4, 'AbsTol', 1e-4,...
%         'InitialStepSize', [],'VariableStepSize', true,'MaxStepSize', 1e1);
% elseif strcmpi(method, 'Merson45')
%     dudt = SSPTools.Steppers.LoadEmbeddedERK('MethodName', 'Merson45',...
%         'ODE', vdp, 'y0', y0, 'RelTol', 1e-4, 'AbsTol', 1e-4,...
%         'InitialStepSize', 0.0001,'VariableStepSize', true,'MaxStepSize', 1e1);
% elseif strcmpi(method, 'HeunEuler')
%     dudt = SSPTools.Steppers.LoadEmbeddedERK('MethodName', 'heuneuler',...
%         'ODE', vdp, 'y0', y0, 'RelTol', 1e-4, 'AbsTol', 1e-4,...
%         'InitialStepSize', [],'VariableStepSize', true,'MaxStepSize', 1e1);
% elseif strcmpi(method, 'Fehlberg12')
%     dudt = SSPTools.Steppers.LoadEmbeddedERK('MethodName', 'Fehlberg12',...
%         'ODE', vdp, 'y0', y0, 'RelTol', 1e-4, 'AbsTol', 1e-4,...
%         'InitialStepSize', [],'VariableStepSize', true,'MaxStepSize', 1e1);
% elseif strcmpi(method, 'BogackiShampine')
%     dudt = SSPTools.Steppers.LoadEmbeddedERK('MethodName', 'BogackiShampine',...
%         'ODE', vdp, 'y0', y0, 'RelTol', 1e-4, 'AbsTol', 1e-4,...
%         'InitialStepSize', [],'VariableStepSize', true,'MaxStepSize', 1e1);
% else
%     dudt = SSPTools.Steppers.EmbeddedERK('A',rk.A, 'b', rk.b,'bhat',rk.bhat', 'p', rk.p, 'phat', rk.p-1, ...
%         'ODE', vdp, 'y0', y0, 'RelTol', 1e-4, 'AbsTol', 1e-4,...
%         'InitialStepSize', 0.001,'VariableStepSize', true,'MaxStepSize', 1e1);
% end
% get the estimated starting step size
%dt_est = dudt.startingStepSize();

dudt = SSPTools.Steppers.LoadEmbeddedERK('MethodName', 'Merson45',...
    'ODE', vdp, 'y0', y0, 'RelTol', 1e-4, 'AbsTol', 1e-4,...
    'InitialStepSize', 0.0001,'VariableStepSize', true,'MaxStepSize', 1e1);
dt = 0.0001;

T = []; DT = []; Y = []; ERR = []; badDT = [];

[t, y, dt, err, bad_dt] = dudt.getState();
T = [T; t]; DT = [DT; dt]; Y = [Y y]; ERR = [ERR; err]; badDT = [badDT; bad_dt];

while t < Tfinal
    
    dudt.takeStep(dt);
    [t, y, nextDt, err, bad_dt] = dudt.getState();
    dt = min(nextDt, Tfinal - t);
    T = [T; t]; DT = [DT; dt]; Y = [Y y]; ERR = [ERR; err]; badDT = [badDT; bad_dt];
        
end


% remove all the rejected Step and Solution

goodSol = isnan(badDT);
badDT = badDT(~goodSol);

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

% % plotting the step size taken at each time step
% subplot(3,1,2)
% title('Accepted Step Size')
% semilogy(T, DT,'-s')
% hold on
% semilogy(T, mean(DT)*ones(size(T)),'--k')
% semilogy(badT, badDT, 'rx')
% ylim([0 5]);
% 
% % plotting the Error at each time-step
% subplot(3,1,3)
% title('Error')
% semilogy(T, ERR,'-.k')
% hold on
% semilogy(T, 1e-4*ones(size(T)),'--')
% ylim([0 5]);
% 
% fprintf(1, '%d accepted + %d rejected\n', dudt.acceptedStep, dudt.rejectedStep);
% %fprintf(1, '%d rejected\n', dudt.rejectedStep);