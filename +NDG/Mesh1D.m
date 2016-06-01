classdef Mesh1D < NDG.Mesh
    
    properties
        
    end
    
    properties (Dependent = true, SetAccess = private)
        
    end
    
    
    methods
        
        function obj = Mesh1D(varargin) % The Constuctor
            % call parent constructor
            obj = obj@NDG.Mesh(varargin{:});
            obj.nNodes = obj.nElements*obj.nSolutionPoints;
            
            
            %function eFaces = get.elementFaces(obj) % Element Faces
            Faces = linspace(obj.domain(1),obj.domain(2),obj.nFaces);
            obj.elementFaces = [Faces(1:end-1); Faces(2:end)];
            
            
            %function eSize = get.elementSize(obj) % Element Size
            %TODO: this is dx
            obj.elementSize = (obj.domain(2)-obj.domain(1))/obj.nElements;
            
            %function eCenter =  get.elementCenter(obj) % Element Center
            obj.elementCenter = (obj.elementFaces(2,:)+obj.elementFaces(1,:))/2;

            %function J = get.Jacobian(obj) % Jacobian dx/dxi
            obj.Jacobian = obj.elementSize/2;
            %end
            
            %function eNodes = get.elementNodes(obj)
            obj.elementNodes = reshape(1:obj.nNodes,obj.nSolutionPoints,obj.nElements);
            %end
            
            %function nodeCoords = get.nodeCoordinates(obj) % x-Grid
            [obj.solPts, obj.weights] = obj.quadFunction(obj.nSolutionPoints);

            [xSPs,xc] = meshgrid(obj.elementCenter,obj.Jacobian*obj.solPts);
            obj.nodeCoordinates = xSPs+xc; % x-Grid            
        end
        
    end % Methods
    
    methods (Access = protected)
        
        
        
    end % protected method
    
end % Class

