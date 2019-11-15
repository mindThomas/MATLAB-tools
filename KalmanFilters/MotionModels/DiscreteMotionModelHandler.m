classdef DiscreteMotionModelHandler
    properties %(SetAccess = private)
        x
    end
	properties (Access = private)
        f
        Q
        Q_sqrt
    end
    methods
        function obj = DiscreteMotionModelHandler(f, x0, Q)                        
            obj.x = x0;
            obj.f = f;
            obj.Q = Q;
            obj.Q_sqrt = chol(Q, 'lower');
        end
        
        function obj = step(obj, vargin)            
            if (nargin == 2)
                u = vargin{1};
            else
                u = [];
            end
            % Generate a realization from the process model based on x and Q
            % Generate noise realization, q[k] ~ N(0, Q)
            q = obj.Q_sqrt * randn([size(obj.Q, 1), 1]);
            % Apply the process model
            obj.x = obj.f(obj.x, u, q); 
        end
        
        function obj = stepDeterministic(obj, vargin)            
            if (nargin == 2)
                u = vargin{1};
            else
                u = [];
            end
            % Zero noise
            q = zeros([size(obj.Q, 1), 1]);
            % Apply the process model
            obj.x = obj.f(obj.x, u, q); 
        end        
        
    end
end