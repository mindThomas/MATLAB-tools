classdef ContinuousMotionModelHandler
    properties %(SetAccess = private)
        x
    end
	properties (Access = private)
        f
        Q
        Q_sqrt
    end
    methods
        function obj = ContinuousMotionModelHandler(f, x0, Q)                        
            obj.x = x0;
            obj.f = f;
            obj.Q = Q;
            obj.Q_sqrt = chol(Q, 'lower');
        end
        
        function obj = step(obj, T, vargin)            
            if (nargin == 3)
                u = vargin{1};
            else
                u = [];
            end
            % For smaller continuous systems, a step from t=t0 to t=t0+T,
            % can be made with Euler discretization of small steps
            steps = 1000;
            ts = T / steps; % do 1.000 steps
            for (i = 1:steps)
                % Generate a realization from the process model based on x and Q
                % Generate noise realization, q[k] ~ N(0, Q)
                q = obj.Q_sqrt * randn([size(obj.Q, 1), 1]);
                % Get the derivative from the motion model
                dx = obj.f(obj.x, u, q); 
                % Compute next step using the Euler method
                %  x[k] = x[k-1] + ts*dx[k-1]
                obj.x = obj.x + ts*dx;
            end
        end
        
        function obj = stepDeterministic(obj, T, vargin)            
            if (nargin == 3)
                u = vargin{1};
            else
                u = [];
            end
            % Zero noise
            q = zeros([size(obj.Q, 1), 1]);
            % For smaller continuous systems, a step from t=t0 to t=t0+T,
            % can be made with Euler discretization of small steps
            steps = 1000;
            ts = T / steps; % do 1.000 steps
            for (i = 1:steps)
                % Get the derivative from the motion model
                dx = obj.f(obj.x, u, q); 
                % Compute next step using the Euler method
                %  x[k] = x[k-1] + ts*dx[k-1]
                obj.x = obj.x + ts*dx;
            end
        end        
        
    end
end