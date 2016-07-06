clear all; close all; clc


% check that the matlab version must be at least R2013b
assert(~verLessThan('matlab', 'R2013b'),'MATLAB R2013b is required');

cellOfTest = {UnitTests.TestOrder; ...  % convergence tests
    UnitTests.TestTVD; ...              % SSP Tests
    UnitTests.TestODEs; ...             % ODE Tets
    UnitTests.TestSSPTools; ...         % run-time Tests
    UnitTests.TestMesh; ...             % Mesh1D Tests
    UnitTests.TestDG;                   % DG Runtime
    UnitTests.TestEmbeddedMethods;       % Embedded-RK
    };

resultTable = [];
for testInd = 1:numel(cellOfTest)
   test = cellOfTest{testInd};
   res = run(test); res.table;
   resultTable = [resultTable;res.table];
end
    
disp(resultTable);