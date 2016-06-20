%Sidafa Conde
%example of Brusselator ODE
% Purpose: recreate result from fig.4.1 in
% Solving Ordinary Differential Equation I
% Hairer. Pg 170

clear all; close all; clc
addpath('../');
addpath('~/Dropbox/SSPTools/');


method = 'DormandPrince54';

t = 0;
y0(1) = 0.994;
y0(2) = 0.0;
y0(3) = 0.0;
y0(4) = -2.00158510637908252240537862224;
Tfinal = 17.0652165601579625588917206249;

aren = TestProblems.ODEs.Aren();

A = zeros(7);
A(2,1)=0.2;
A(3,1)=3.0/40.0;
A(3,2)=9.0/40.0;
A(4,1)=44.0/45.0;
A(4,2)=-56.0/15.0;
A(4,3)=32.0/9.0;
A(5,1)=19372.0/6561.0;
A(5,2)=-25360.0/2187.0;
A(5,3)=64448.0/6561.0;
A(5,4)=-212.0/729.0;
A(6,1)=9017.0/3168.0;
A(6,2)=-355.0/33.0;
A(6,3)=46732.0/5247.0;
A(6,4)=49.0/176.0;
A(6,5)=-5103.0/18656.0;
A(7,1)=35.0/384.0;
A(7,3)=500.0/1113.0;
A(7,4)=125.0/192.0;
A(7,5)=-2187.0/6784.0;
A(7,6)=11.0/84.0;

b = zeros(7,1);
b(1) = 35/384;
b(3) = 500/1113;
b(4) = 125/192;
b(5) = -2187/6784;
b(6) = 11/84;

bhat = zeros(7,1);
bhat(1) = 5179/57600;
bhat(3) = 7571/16695;
bhat(4) = 393/640;
bhat(5) = -92097/339200;
bhat(6) = 187/2100;
bhat(7) = 1/40;

c = sum(A,2);

E = zeros(size(c));
E(1)=71.0/57600.0;
E(3)=-71.0/16695.0;
E(4)=71.0/1920.0;
E(5)=-17253.0/339200.0;
E(6)=22.0/525.0;
E(7)=-1.0/40.0;

D = zeros(size(E));
D(1)=-12715105075.0/11282082432.0;
D(3)=87487479700.0/32700410799.0;
D(4)=-10690763975.0/1880347072.0;
D(5)=701980252875.0/199316789632.0;
D(6)=-1453857185.0/822651844.0;
D(7)=69997945.0/29380423.0;

p = 4; phat = 5;

dudt = SSPTools.Steppers.EmbeddedERK('A',A,'b', b, 'bhat',bhat, 'p',p, 'phat', phat,...
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