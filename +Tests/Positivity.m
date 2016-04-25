classdef Positivity < Tests.Test
   
    properties
        name;
    
    end
    
    methods
        function obj = Positivity(varargin)  %constructor
            obj = obj@Tests.Test(varargin{:});
        end
        
        
        function run(varargin) end
        
        function [ output ] = run_test(varargin) end
    end
    
    methods (Acsess = protected)
        function log(obj, varargin) end 
    end
    
end