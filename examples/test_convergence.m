clear all; close all; clc

addpath('../')
N = 16;
Tfinal = 0.4;

%testing = 'ERK';
%testing = 'DIRK';
testing = 'IMEXRK';
%testing = 'IMEXRK-SSP11';


y0 = @(x) sin(x);

% LNL-IMEX RK
pex = 2;
pim = 2;
plin = 2;
s = 5;
k = 0;

file = sprintf('~/Dropbox/imex-linear/src/butcher-optimization/Method/DIRK/G/Pex%d/Pim%d/Plin%d/S%d/K',...
    pex, pim, plin,s);
file = [file num2str(k) '/'];
files = dir([file '*.mat']);
method = files(end).name;
rk = load([file method]);
A = rk.A; b = rk.b; s = rk.s;
At = rk.At; bt = rk.bt;

A = [0]; b = [1]; s = 1; %TODO: infert s from size(A,1)
At = [1]; bt = [1];

exp_pro = TestProblems.PDEs.LinearAdvection('a', 1);
exp_pro = TestProblems.PDEs.Burgers();
imp_pro = TestProblems.PDEs.BuckleyLeverett();

dfdx = SSPTools.Discretizers.Spectral('derivativeOrder',1, 'N', N);

if strcmpi(testing, 'erk')
    dudt = SSPTools.Steppers.ERK('A', A, 'b',b, 's', s,...
        'dfdx', dfdx, 'ExplicitProblem', imp_pro, 'y0', y0);
elseif strcmpi(testing,'dirk')
    dudt = SSPTools.Steppers.DIRK('A', At, 'b',bt, 's', s,...
        'dfdx', dfdx, 'ExplicitProblem', imp_pro, 'y0', y0);
elseif strcmpi(testing, 'imexrk')
    dudt = SSPTools.Steppers.IMEXRK('A', A, 'b',b, 's', s, 'At', At, 'bt', bt,...
        'dfdx', dfdx, 'ExplicitProblem', exp_pro, 'ImplicitProblem', imp_pro,...
        'dgdx', dfdx, 'y0',y0);
elseif strcmpi(testing, 'imexrk-ssp11')
    dudt = SSPTools.Steppers.LoadIMEX('MethodName', 'Stormer-Verlet','dfdx', dfdx, 'ExplicitProblem', exp_pro, 'ImplicitProblem', imp_pro,...
        'dgdx', dfdx, 'y0',y0);
end

convergencePDE = Tests.Convergence('integrator', dudt,'Tfinal', Tfinal,...
    'CFL', (1/2).^(1:5));

convergencePDE.run();
convergencePDE.complete();
%convergencePDE.getOrder('L2')
%convergencePDE.getOrder('L1')
%convergencePDE.getOrder('Linf')