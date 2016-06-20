classdef EmbeddedERK < SSPTools.Steppers.EmbeddedRK
    
    properties
        %bhat; % embedding weight vector
    end
    
    properties ( Access = private)
        variableStep = true;    
    end
    
    methods
        
        function obj = EmbeddedERK(varargin)
            obj = obj@SSPTools.Steppers.EmbeddedRK(varargin{:});
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'name','Embedded-ERK');
            
            p.parse(varargin{:});
            
            %obj.variableStep = p.Results.VariableStepSize;
            
            obj.name = p.Results.name;

        end % end constructor
        
        
    end
    
    methods %( Access = protected )
        
        
        function [y, dt] = takeStep(obj, dt)
            % function [y, dt] = takeStep(dt)
            % returns the new solution (y) and time-step taken (dt)
            
            u0 = obj.u0;
            obj.Y(:,1) = u0;
            obj.dt_ = dt;
            obj.nstep = obj.nstep + 1;
            
            % intermediate stage value
            for i = 2:obj.s
                temp = u0;
                for j = 1:i-1
                    temp = temp + dt*obj.A(i,j)*obj.L(dt + obj.c(j), obj.Y(:,j));
                end
                obj.Y(:,i) = temp;
            end
            
            obj.nfcn = obj.nfcn + (obj.s-1);
            
            y = u0;
            for i = 1:obj.s
                y = y + dt*obj.b(i)*obj.L(dt + obj.c(i), obj.Y(:,i));
            end
            
            % compute the embedding solution
            if obj.isEmbedded % just to make sure
                yhat = u0;
                for i = 1:obj.s
                    yhat = yhat + dt*obj.bhat(i)*obj.L(dt + obj.c(i), obj.Y(:,i));
                end
                
                [y, t] = obj.stepSizeControl(dt, y, yhat);
            else
                obj.nextDt = dt;
                t = obj.t + dt;
            end
            
            obj.u0 = y;
            obj.t  = t;
        end

   
    end
    
    methods ( Access = protected )
        
        

    end
    
    
end
