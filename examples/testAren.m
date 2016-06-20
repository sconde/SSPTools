%Sidafa Conde
%example of Brusselator ODE
% Purpose: recreate result from fig.4.1 in
% Solving Ordinary Differential Equation I
% Hairer. Pg 170

clear all; close all; clc
addpath('../');
addpath('~/Dropbox/SSPTools/');


%method = 'DormandPrince54';
%method = 'Merson45';
%method = 'Zonneveld43';
%method = 'Felhberg45';
%method = '38rule43';

t = 0;
y0(1) = 0.994;
y0(2) = 0.0;
y0(3) = 0.0;
y0(4) = -2.00158510637908252240537862224;
Tfinal = 17.0652165601579625588917206249;

aren = TestProblems.ODEs.Aren();

dudt = SSPTools.Steppers.LoadEmbeddedERK('MethodName',method,...
    'ODE', aren, 'y0', y0, 'RelativeTolerance', 1e-7, 'AbsoluteTolerance', 1e-7,...
    'InitialStepSize', [],'VariableStepSize', true, 'Tfinal', Tfinal);

T = []; DT = []; Y = []; ERR = []; badDT = [];

[t, y, dt, err, bad_dt] = dudt.getState();
T = [T; t]; DT = [DT; dt]; Y = [Y y]; ERR = [ERR; err]; badDT = [badDT; bad_dt];

while t < Tfinal
    
    dudt.takeStep(dt);
    [t, y, nextDt, err, bad_dt] = dudt.getState();
    
    dt = min(nextDt, Tfinal - t);
    T = [T; t]; DT = [DT; dt]; Y = [Y y]; ERR = [ERR; err]; badDT = [badDT; bad_dt];
end

% fprintf(1, 'rtol= %12.5e, fcn = %d, step = %d, accpt = %d, reject = %d\n',...
%     dudt.absTol, dudt.nfcn, dudt.nstep, dudt.naccpt, dudt.nrejct);
dudt.summary();

% T = T(2:end); DT = DT(2:end); ERR = ERR(2:end); badDT = badDT(2:end);
% % remove all the rejected Step and Solution
%
% goodSol = isnan(badDT);
% badDT = badDT(~goodSol);
%
%
% [ismem, ind] = ismember(badDT,DT);
% DT(ind) = [];
% ERR(ind) = [];
% Y = Y';
% Y = Y(goodSol,:);
% badT = T(~goodSol);
% T = T(goodSol);
%
% figure('Position',[100 100 1200 400]);
% % plotting the solution
% subplot(3,1,1);
% title('Solution')
% plot(T, Y(:,2), '-or')
% hold on
% plot(T, Y(:,1), '-pk')
% xlim([0 Tfinal]);
%
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
% fprintf(1, '%d Function-Evaluation, %d Steps( %d accepted + %d rejected)\n',...
%     dudt.nfcn, dudt.nstep, dudt.naccpt, dudt.nrejct);
% %fprintf(1, '%d rejected\n', dudt.rejectedStep);