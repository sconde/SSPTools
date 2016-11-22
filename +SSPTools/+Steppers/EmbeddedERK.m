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
            t_ = obj.t;
            obj.Y(:,1) = u0;
            obj.dt_ = dt;
            obj.Fvec(:,1) = obj.L(t_ , u0);
            obj.nstep = obj.nstep + 1;
            
            % intermediate stage value
            for i = 2:obj.s
                temp = u0;
                for j = 1:i-1
                    temp = temp + dt*obj.A(i,j)*obj.Fvec(:, j);
                end
                obj.Fvec(:, i) = obj.L(t_ + obj.c(i)*dt, temp);
            end
            
            obj.nfcn = obj.nfcn + (obj.s-1);
            
            % combine
            y = u0 + dt*obj.Fvec*obj.b(:);
            
            % compute the embedding solution
            if obj.isEmbedded % just to make sure
                % combine
                yhat = u0 + dt*obj.Fvec*obj.bhat(:);
                
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
