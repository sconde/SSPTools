function rk = loadMethod(pex, pim, plin, s, k)


file = sprintf('~/Dropbox/imex-linear/src/butcher-optimization/Method/DIRK/G/Pex%d/Pim%d/Plin%d/S%d/K',...
    pex, pim, plin,s);
file = [file num2str(k) '/'];
files = dir([file '*.mat']);
method = files(end).name;
rk = load([file method]);
end