classdef Test < handle
   
    properties
        name = 'Convergence';
    end
    
    methods
        function obj = Test(varargin) end
        
        function run(varargin) end
        
        function [ output ] = run_test(varargin) end
    end
    
    methods %(Acsess = protected)
        function log(obj, varargin) end 
    end
    
end
