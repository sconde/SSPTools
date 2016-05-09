clear all; close all; clc

addpath('../');
N = 300;

%testing = 'ERK';
%testing = 'DIRK';
testing = 'IMEXRK';
%testing = 'IMEXRK-SSP11';


y0 = @(x) heaviside(x - (ceil((x+1)/2) -1)*2);

% LNL-IMEX RK
pex = 3;
pim = 2;
plin = 4;
s = 5;
k = 1;

file = sprintf('~/Dropbox/imex-linear/src/butcher-optimization/Method/DIRK/G/Pex%d/Pim%d/Plin%d/S%d/K',...
    pex, pim, plin,s);
file = [file num2str(k) '/'];
files = dir([file '*.mat']);
method = files(end).name;
rk = load([file method]);
A = rk.A; b = rk.b; s = rk.s;
At = rk.At; bt = rk.bt;
%[rk.r rk.rt]


imp_pro = TestProblems.PDEs.LinearAdvection('a', 1);
exp_pro = TestProblems.PDEs.Burgers();
%imp_pro = TestProblems.PDEs.BuckleyLeverett();

dfdx = SSPTools.Discretizers.FiniteDifference('N', N, 'domain', [-1, 1],'bc','periodic');

if strcmpi(testing, 'erk')
    dudt = SSPTools.Steppers.ERK('A', rk.A, 'b',rk.b, 'p', rk.pex,'plin', rk.plin,...
       'dfdx', dfdx, 'ExplicitProblem', imp_pro, 'y0', y0);
    
%     dudt = SSPTools.Steppers.LoadERK('MethodName', 'FE',...
%         'dfdx', dfdx, 'ExplicitProblem', imp_pro, 'y0', y0);
elseif strcmpi(testing,'dirk')

    dudt = SSPTools.Steppers.DIRK('A', At, 'b',bt, 'p', rk.pim,'plin', rk.plin,...
        'dfdx', dfdx, 'ExplicitProblem', imp_pro, 'y0', y0);
elseif strcmpi(testing, 'imexrk')
    dudt = SSPTools.Steppers.IMEXRK('A', A, 'b',b, 'At', At, 'bt', bt,...
        'p', rk.pex, 'pim',rk.pim, 'plin',rk.plin,...
        'dfdx', dfdx, 'ExplicitProblem', exp_pro, 'ImplicitProblem', imp_pro,...
        'dgdx', dfdx, 'y0',y0);
elseif strcmpi(testing, 'imexrk-ssp11')
    dudt = SSPTools.Steppers.LoadIMEX('MethodName', 'IMEX1','dfdx', dfdx, 'ExplicitProblem', exp_pro, 'ImplicitProblem', imp_pro,...
        'dgdx', dfdx, 'y0',y0);
end

tvbPDE = Tests.SSP('integrator', dudt,'TVB',true,'CFLRefinement',0.1,...
    'CFLMAX',2,'CFL',0.8);
tvdPDE = Tests.SSP('integrator', dudt,'TVD',true,'CFLRefinement',0.001,...
    'CFLMAX',1.1,'CFL',0.95);

keyboard

% tvbPDE.run();
% tvbPDE.plotSolution()
% 
% keyboard
tvdPDE.run();
tvdPDE.plotSolution()