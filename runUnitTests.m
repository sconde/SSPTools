clear all; close all; clc

cellOfTest = {UnitTests.TestOrder; ...  % convergence tests
    UnitTests.TestTVD; ...              % SSP Tests
    UnitTests.TestODEs; ...             % ODE Tets
    UnitTests.TestSSPTools; ...         % run-time Tests
    UnitTests.TestMesh; ...             % Mesh1D Tests
    UnitTests.TestDG;                   % DG Runtime
    };

for testInd = 1:numel(cellOfTest)
   test = cellOfTest{testInd};
   res = run(test); res.table
end
    
% % run all the convergence tests
% test1 = UnitTests.TestOrder;
% res = run(test1); res.table
% 
% % run all the ssp tests
% test2 = UnitTests.TestTVD;
% res2 = run(test2); res2.table
% 
% % run all the ODE tests
% test3 = UnitTests.TestODEs;
% res3 = run(test3); res3.table
% 
% % run all run-time test 
% test4 = UnitTests.TestSSPTools;
% res4 = run(test4); res4.table
% 
% % run Mesh1D tests
% test5 = UnitTests.TestMesh;
% res5 = run(test5); res5.table
% 
% % run TestDG with FE
% test6 = UnitTests.TestDG;
% res6 = run(test6); res6.table