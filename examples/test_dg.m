clear all; close all; clc

addpath('../');
addpath('~/Documents/NDG/');

%testing = 'ERK';
testing = 'DIRK';
Tfinal = 1;
dt = 0.01;

%y0 = @(x) heaviside(x - (ceil((x+1)/2) -1)*2);
y0 = @(x) sin(pi*x);

nE = 15;    % number of elements
K  = 5;     % degree of accuracy

exp_pro = TestProblems.PDEs.LinearAdvection('a',1);

xmesh = NDG.Mesh1D('Domain', [-1 1], 'NumberElements', nE,...
    'SolutionDegree', K, 'QuadratureType', 'LGL');

% ndg_mesh = mesh1d([0 2*pi],nE,'LGL',K);
% 
% ndg_dx = ndg_mesh.elementSize;     dx = xmesh.elementSize;      assert(isequaln(ndg_dx,dx));
% ndg_J = ndg_mesh.Jacobian;         J = xmesh.Jacobian;          assert(isequaln(ndg_J,J));
% ndg_x = ndg_mesh.nodeCoordinates;  x = xmesh.nodeCoordinates;   assert(isequaln(ndg_x, x));
% ndg_w = ndg_mesh.weights';         w = xmesh.weights';
% ndg_xc = ndg_mesh.elementCenter;   xc = xmesh.elementCenter;


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