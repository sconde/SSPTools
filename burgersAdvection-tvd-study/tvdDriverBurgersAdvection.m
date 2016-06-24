clear all; close all; clc;

global dudt 

rk = loadMethod(2,2,2,2,1);
r = rk.r;
K = rk.k;

a = 1;
N = 601;

y0 = @(x) heaviside(x - (ceil((x+1)/2) -1)*2);


addpath('~/Dropbox/SSPTools/')

domain = [-1 1];

% define the advection flux, a = wave speed
imp_pro = TestProblems.PDEs.LinearAdvection('a',a);

% define the burgers flux
exp_pro = TestProblems.PDEs.Burgers();

% define the semi-discrete explcicit problem
dfdx = SSPTools.Discretizers.FiniteDifference('derivativeOrder',1, 'N', N,...
    'Problem', exp_pro, 'domain', domain);

% define the semi-discrete implicit problem
dgdx = SSPTools.Discretizers.FiniteDifference('derivativeOrder',1, 'N', N,...
    'Problem', imp_pro, 'domain', domain);

dudt = SSPTools.Steppers.IMEXRK('A', rk.A, 'b',rk.b, 'At', rk.At, 'bt', rk.bt,...
    'dfdx', dfdx,'dgdx', dgdx, 'y0', y0);


startingr = min(r,.005);
sp = min(r/10,.005);

lambda = sort(startingr + rand(1,20)*(r + 0.5 - startingr));
lambda = linspace(startingr, r + 0.5, 10);

cfl_refinement = max(diff(lambda))

L = [];
VV = [];
V = 0;
coarse = 1;
VV_ = [];


while cfl_refinement > 10e-10
    Violation_ = burgersAdvection( lambda);
    
    L = [L,lambda];
    VV_ = [VV_, Violation_];
    
    % find the first violation
    ind = find(log10(Violation_) > -12, 1, 'first');
    
    if isempty(ind);  %If the observed CFL is outside original range
        maxL = max(lambda);
        ind = find(lambda == maxL, 1, 'first');
        lambda = linspace(maxL-cfl_refinement, maxL + 2*cfl_refinement,10);
    else
        Ltemp = sort(L);
        ind_ = find(Ltemp == lambda(ind),1,'first');
        newL = Ltemp(ind_);
        lambda = linspace(newL-2*cfl_refinement,newL + 2*cfl_refinement,10);
    end
    cfl_refinement = max(diff(lambda))
    
end

log10VV = log10(VV_);

% chop off the extreme values
extremeInd = log10VV >= 0;
log10VV(extremeInd) = 0;

% get the index of last stable cfl (this is the observed SSP)
indGood = log10VV <= -14;
goodSSPInd = find(indGood, 1, 'last');
observedSSP = L(goodSSPInd);

plot(L, log10VV, 'kx', 'markersize', 8);
hold on
plot(L(goodSSPInd), log10VV(goodSSPInd),  'pb', 'markersize', 20)
title(sprintf('r = %5.4f, observed = %5.4f', r, observedSSP));
xlabel('$\frac{\Delta t}{\Delta x}$', 'Interpreter' ,'latex','FontSize',20)
ylabel('Log10 of Maximum Violation of Total Variation','FontSize',14)
%ylabel('log_{10} \max{\| u^{n+1} \|_{TV} - \| u^{n}\|_{TV}}','FontSize',14)
