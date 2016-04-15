classdef RK < handle
    % Integrator.m
    %
    % This is the base class for all numerical time-stepping methods in the SSP_Tools
    % package. It provides a simple external interface for communicating with a time-
    % stepping method along with basic support mechanisms for evaluating the
    %
    %
    %
    %
    properties
        name;			% Name of time-stepping method.
        dydt;
        y0;
    end
    
    methods
        function obj = RK(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParamValue(p,'name','MSRK');
            p.parse(varargin{:});
            obj.name = p.Results.name;
        end
    end
    
    methods %( Access = protected )
        function [y] = takeStep(obj, dt)
        
        end
        
    end
    
end