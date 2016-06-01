clear all; close all; clc

addpath('../');
addpath('~/Documents/NDG/');

%testing = 'ERK';
testing = 'DIRK';
Tfinal = 2;
dt = 0.01;

y0 = @(x) heaviside(x - (ceil((x+1)/2) -1)*2);
y0 = @(x) sin(pi*x);

nE = 20;    % number of elements
K  = 2;     % degree of accuracy

%exp_pro = TestProblems.PDEs.LinearAdvection('a',1);
exp_pro = TestProblems.PDEs.Burgers();

xmesh = NDG.Mesh1D('Domain', [-1 1], 'NumberElements', nE,...
    'SolutionDegree', K, 'QuadratureType', 'LGL');


dg = NDG.NDG('Mesh', xmesh, 'Problem', exp_pro);

dudt = SSPTools.Steppers.LoadERK('MethodName', 'FE',...
    'dfdx', dg, 'y0', y0);


line1 = plot(dg.x(:), dudt.y0(dg.x(:)),'-r','linewidth',2);
axis([-1 1 -1 1]);

t = 0;
while t < Tfinal
    
    ynew = dudt.takeStep(dt);
    t = t+ dt
    
    set(line1, 'ydata', ynew);
    drawnow;
    pause(0.1);
    
end