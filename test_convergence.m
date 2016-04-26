clear all; close all; clc


N = 8;
Tfinal = 1;

testing = 'ERK';
%testing = 'DIRK';


y0 = @(x) sin(x);
f = @(t, u) u;

% LNL-IMEX RK
pex = 3;
pim = 3;
plin = 4;
s = 6;
k = 1;

file = sprintf('~/Dropbox/imex-linear/src/butcher-optimization/Method/DIRK/G/Pex%d/Pim%d/Plin%d/S%d/K',...
    pex, pim, plin,s);
file = [file num2str(k) '/'];
files = dir([file '*.mat']);
method = files(end).name;
rk = load([file method]);
A = rk.A; b = rk.b; s = rk.s;
At = rk.At; bt = rk.bt;

% A = [0 0;1 0]; b = [1 0]; s = 2; %TODO: infert s from size(A,1)
% At = [0 0;0 1]; bt = [0 1];

imp_pro = TestProblems.PDEs.LinearAdvection('a', 1);

dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, 'N', N);

if strcmpi(testing, 'erk')
    dudt = SSPTools.Steppers.ERK('A', A, 'b',b, 's', s,...
        'dfdx', dfdx, 'ExplicitProblem', imp_pro, 'y0', y0);
else
    dudt = SSPTools.Steppers.DIRK('A', At, 'b',bt, 's', s,...
        'dfdx', dfdx, 'ExplicitProblem', imp_pro, 'y0', y0);
end

convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal', Tfinal,...
    'CFL', (1/2).^(1:5));

convergencePDE.run();
convergencePDE.complete();