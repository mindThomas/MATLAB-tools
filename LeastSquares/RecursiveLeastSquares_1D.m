classdef RecursiveLeastSquares_1D
    % Different implementations of ordinary recursive least squares
    % min ||a*x - y||^2
    % where || is the 2-norm (Euclidean distance)
    % with an recursively incoming value of a_n and measurement y_n
    properties
        x_est
    end
    properties (Access = private)
        a_squared
        a_squared_inv
        ay
    end    
    
    methods
        function obj = RecursiveLeastSquares_1D(a_init, y_init)            
            % Pseudo inverse is given as
            % Ainv = inv(A'*A) * A';
            % x = Ainv * y;
            obj.a_squared = a_init^2;
            obj.ay = a_init * y_init;
            obj.a_squared_inv = 1/obj.a_squared;
            obj.x_est = 0;
        end   
        
        function obj = AddMeasurement(obj, a_n, y_n)
            % a_n should be the new value of 'a' to append
            % http://pfister.ee.duke.edu/courses/ece586/ex_proj_2008.pdf            
            % Use normal linear algebra to extend A'*A and A'*y
            obj.a_squared = obj.a_squared + a_n^2;
            obj.ay = obj.ay + a_n*y_n;
        end
        
        function obj = AddMeasurementAndComputeEstimate(obj, a_n, y_n)
            % a_n should be the new value of 'a' to append
            % http://pfister.ee.duke.edu/courses/ece586/ex_proj_2008.pdf            
            % Use normal linear algebra to extend A'*A and A'*y
            obj.a_squared = obj.a_squared + a_n^2;
            obj.ay = obj.ay + a_n*y_n;
            
            % Compute estimate
            obj.x_est = obj.ay / obj.a_squared;
        end
        
        function x_est = ComputeEstimate(obj)            
            % Compute estimate
            x_est = obj.ay / obj.a_squared;           
        end
    end
end
