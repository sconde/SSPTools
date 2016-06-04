clear all; close all; clc

cellOfTest = {UnitTests.TestOrder; ...  % convergence tests
    UnitTests.TestTVD; ...              % SSP Tests
    UnitTests.TestODEs; ...             % ODE Tets
    UnitTests.TestSSPTools; ...         % run-time Tests
    UnitTests.TestMesh; ...             % Mesh1D Tests
    UnitTests.TestDG;                   % DG Runtime
    };

resultTable = [];
for testInd = 1:numel(cellOfTest)
   test = cellOfTest{testInd};
   res = run(test); res.table;
   resultTable = [resultTable;res.table];
end
    
disp(resultTable);