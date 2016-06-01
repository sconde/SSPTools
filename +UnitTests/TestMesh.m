classdef TestMesh < matlab.unittest.TestCase
    
    properties
        
    end
    
    % Test Method Block
    methods (Test)
        % includes unit test functions
        
        function testMesh1DLGL(testCase)
            % test for completition
            
            import matlab.unittest.TestCase
            
            nE = 15;    % number of elements
            K  = 2;     % degree of accuracy
            
            xmesh = NDG.Mesh1D('Domain', [-1 1], 'NumberElements', nE,...
                'SolutionDegree', K, 'QuadratureType', 'LGL');
            
            testCase.assertTrue(true)
            assertTrue(testCase, true)
            
        end
        
               function testMesh1DLegendre(testCase)
            % test for completition
            
            import matlab.unittest.TestCase
            
            nE = 15;    % number of elements
            K  = 2;     % degree of accuracy
            
            xmesh = NDG.Mesh1D('Domain', [-1 1], 'NumberElements', nE,...
                'SolutionDegree', K, 'QuadratureType', 'Legendre');
            
            testCase.assertTrue(true)
            assertTrue(testCase, true)
            
        end
        
        
        
    end
end
