clear all; close all; clc

% run all the convergence tests
test1 = UnitTests.TestOrder;
res = run(test1); res.table

% run all the ssp tests
test2 = UnitTests.TestTVD;
res2 = run(test2); res2.table

% run all the ODE tests
test3 = UnitTests.TestODEs;
res3 = run(test3); res3.table