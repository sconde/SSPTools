classdef NDG < NDG.NDGBase
    
    properties
        
    end
    
    properties (Access = protected)
        Lift;
        Dr;
    end
    
    methods
        function obj = NDG(varargin)
            % TODO: determine the default values?
            obj = obj@NDG.NDGBase(varargin{:});
            
            p = inputParser;
            p.KeepUnmatched = true;
            
            % build the Lift Operator
            Emat = zeros(obj.mesh.solutionDegree + 1, 2);
            Emat(1,1) = 1; Emat(obj.mesh.solutionDegree + 1, 2) = 1;
            obj.Lift = obj.Vandermonde2*(obj.Vandermonde2'*Emat);
            obj.Dr = obj.nodalCoefDiffMatrix;
            obj.dx = obj.mesh.elementSize;
            obj.x = obj.mesh.nodeCoordinates;
        end % NDG constructor
        
        function [y] = L(obj, t, u)
            
            u = reshape(u, size(obj.mesh.nodeCoordinates));
            % first call the residual function
            dF = obj.residual(t, u, obj.problem.f, obj.problem.em, obj.Lift, obj.Dr);
            y = -dF/obj.mesh.Jacobian;
            y = y(:)*obj.dx; % make a column vector
        end
    end
    
end
