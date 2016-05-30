clear all; close all; clc

addpath('../');
N = 300;

testing = 'ERK';
%testing = 'DIRK';

y0 = @(x) heaviside(x - (ceil((x+1)/2) -1)*2);

imp_pro = TestProblems.PDEs.LinearAdvection('a', 1);

dfdx = SSPTools.Discretizers.FiniteDifference('N', N, 'domain', [-1, 1],...
    'bc','periodic','Problem', imp_pro);

if strcmpi(testing, 'erk')
    
    
    dudt = SSPTools.Steppers.LoadERK('MethodName', 'FE',...
        'dfdx', dfdx, 'y0', y0);
    
    tvdPDE = Tests.SSP('integrator', dudt,'TVD',true,'CFLRefinement',0.01,...
        'CFLMAX',1.1,'CFL',0.85);
    
elseif strcmpi(testing,'dirk')
    
    dudt = SSPTools.Steppers.LoadDIRK('MethodName', 'MidPoint',...
        'dfdx', dfdx, 'y0', y0);
    
    tvdPDE = Tests.SSP('integrator', dudt,'TVD',true,'CFLRefinement',0.1,...
        'CFLMAX',2,'CFL',0.85);
end

tvbPDE = Tests.SSP('integrator', dudt,'TVB',true,'CFLRefinement',0.1,...
    'CFLMAX',2,'CFL',0.8);

tvdPDE.run();
tvdPDE.plotSolution();
