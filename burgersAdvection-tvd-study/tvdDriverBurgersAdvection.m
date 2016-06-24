clear all; close all; clc;

rk = loadMethod(2,2,2,2,10);
r = rk.r;
K = rk.k;


A = [1 5 10 15];
CFL = cell(1, numel(A));
TV = cell(1, numel(A));

MARKERCell = {'pk','>c','^r','sb','<m','hr','pk',...
    '^c','^r', 'sb', '<m'};
fig = figure('position',[100 100 850 600]); clf;
LISTOFLEG = cell(size(A));
for i = 1:numel(A)
    [  cfl, tv , ssp] =  burgersAdvection( rk, A(i) ) ;
    CFL{i} = cfl;
    TV{i} = tv;
    plot(cfl, tv, MARKERCell{i}); hold on
    LISTOFLEG{i} = sprintf('a = %d, ssp = %5.4f', A(i), ssp);
end

hleg = legend(LISTOFLEG, 'Location','NorthEast');
set(hleg, 'FontSize', 12);

title(sprintf('k = %3.2f, r = %5.4f', rk.k, rk.r));