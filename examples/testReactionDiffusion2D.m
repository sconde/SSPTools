clear all; close all; clc

addpath('../')
A = [0 0;1 0]; b = [1 0]; s = 2; %TODO: infert s from size(A,1)
At = [0 0;0 1]; bt = [0 1];
N = 16;

dt = 0.01;
Tfinal = 1;
t = 0;

testing = 'ERK';

% define the initial conditions

a = 1;
b = 1;
c = 1;
d = 1;

y0 = { @(x, y, t) sin(pi.*x).*sin(pi.*y).*cos(pi.*t);
    @(x, y, t) sin(pi.*x).*sin(pi.*y).*cos(pi.*t/2)};

forcingFnc = {@(t,x,y,a, b, c, d) -pi.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y)+a.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y).*1.973920880217872e1-c.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y);
    @(t,x,y,a, b, c, d)pi.*sin(pi.*t.*(1.0./2.0)).*sin(pi.*x).*sin(pi.*y).*(-1.0./2.0)+b.*cos(pi.*t.*(1.0./2.0)).*sin(pi.*x).*sin(pi.*y).*1.973920880217872e1-d.*cos(pi.*t.*(1.0./2.0)).*sin(pi.*x).*sin(pi.*y)};
    

imp_pro = TestProblems.PDEs.ReactionDiffusion2D('a', 1, 'b', 0.5, 'c', 1,...
    'd', -1,'SourceTerm', forcingFnc);


dfdx = SSPTools.Discretizers.FiniteDifference('derivativeOrder',2, 'N', N,...
    'Problem', imp_pro, 'Domain', [-1 1], 'bc', 'periodic', 'OrderAccuracy', 4, 'Direction', 'CD', 'Dimension', 2);

dudt = SSPTools.Steppers.LoadERK('MethodName','FE',...
    'dfdx', dfdx,'y0', y0);


% % need to initialize DIRK stepper better
% dudt = SSPTools.Steppers.LoadDIRK('MethodName','SDIRK22',...
%     'dfdx', dfdx,'y0', y0);


[t, y] = dudt.getState();

u1 = y(:,1); u2 = y(:,2);
U1 = reshape(u1, size(dfdx.x));
U2 = reshape(u2, size(dfdx.x));

subplot(2,1,1);
ln1 = mesh(dfdx.x, dfdx.y, U1);

subplot(2,1,2);
ln2 = mesh(dfdx.x, dfdx.y, U2);

%line1 = plot(dfdx.x, dudt.y0(dfdx.x),'-r','linewidth',2);
%axis([-1 1 0 1]);


while t < Tfinal
    
    dudt.takeStep(dt);
    
    [t, ynew] = dudt.getState();
    t = t+ dt
        
    u1 = ynew(:,1); u2 = ynew(:,2);
    U1 = reshape(u1, size(dfdx.x));
    U2 = reshape(u2, size(dfdx.x));
    
    set(ln1, 'zdata', U1);
    set(ln2, 'zdata', U2);
    drawnow;
    pause(0.1);
    
end