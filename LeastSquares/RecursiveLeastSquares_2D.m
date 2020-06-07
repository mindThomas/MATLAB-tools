classdef RecursiveLeastSquares_2D
    % Different implementations of ordinary recursive least squares
    % min ||y - A*x||^2    
    % where || is the 2-norm (Euclidean distance)
    % with a recursively incoming 2-element row of A_n and measurement y_n
    properties
        x_est
    end
    properties (Access = private)
        AtA
        AtA_inv
        AtY
    end    
    
    methods
        function obj = RecursiveLeastSquares_2D(A_init, y_init)  
            if (size(y_init,1) ~= 2 ||size(y_init,2) ~= 1)
                error('Incorrect size of y_init. Needs to be a column vector with 2 elements.');
            end
            if (size(A_init,1) ~= 2 ||size(A_init,2) ~= 2)
                error('Incorrect size of A_init. Needs to be a 2x2 matrix.');
            end            
            % Pseudo inverse is given as
            % Ainv = inv(A'*A) * A';
            % x = Ainv * y;
            
            % Compute A'*A
            obj.AtA = zeros(2,2);
            obj.AtA(1,1) = A_init(1,1)*A_init(1,1) + A_init(2,1)*A_init(2,1);
            obj.AtA(1,2) = A_init(1,1)*A_init(1,2) + A_init(2,1)*A_init(2,2);
            obj.AtA(2,1) = A_init(1,2)*A_init(1,1) + A_init(2,2)*A_init(2,1);
            obj.AtA(2,2) = A_init(1,2)*A_init(1,2) + A_init(2,2)*A_init(2,2);
            
            % Compute A'*y
            % obj.AtY = A_init' * y_init;
            obj.AtY = zeros(2,1);
            obj.AtY(1) = A_init(1,1)*y_init(1) + A_init(2,1)*y_init(2);
            obj.AtY(2) = A_init(1,2)*y_init(1) + A_init(2,2)*y_init(2);
            
            % Compute inv(A'*A)
            det = 1 / (obj.AtA(1,1)*obj.AtA(2,2) - obj.AtA(1,2)*obj.AtA(2,1));
            obj.AtA_inv = zeros(2,2);
            obj.AtA_inv(1,1) = det * obj.AtA(2,2);
            obj.AtA_inv(1,2) = -det * obj.AtA(1,2);
            obj.AtA_inv(2,1) = -det * obj.AtA(2,1);
            obj.AtA_inv(2,2) = det * obj.AtA(1,1);
            
            % Compute estimate, AtA_inv * obj.AtY
            obj.x_est = zeros(size(A_init,2), 1);
            obj.x_est(1) = obj.AtA_inv(1,1)*obj.AtY(1) + obj.AtA_inv(1,2)*obj.AtY(2);
            obj.x_est(2) = obj.AtA_inv(2,1)*obj.AtY(1) + obj.AtA_inv(2,2)*obj.AtY(2);                        
        end   
        
        function obj = AddMeasurement(obj, A_n, y_n)
            % A_n should be the new row of A to append
            if (size(A_n,1) ~= 1)
                error('A_n should be the new row of A to append');
            end
            if (size(A_n,2) ~= length(obj.AtY))
                error('Size of A_n does not match number of parameters to estimate');
            end
            % http://pfister.ee.duke.edu/courses/ece586/ex_proj_2008.pdf            
            % Use normal linear algebra to extend A'*A and A'*y
            % A'*A += A_n'*A_n;
            obj.AtA(1,1) = obj.AtA(1,1) + A_n(1)*A_n(1);
            obj.AtA(1,2) = obj.AtA(1,2) + A_n(1)*A_n(2);
            obj.AtA(2,1) = obj.AtA(2,1) + A_n(2)*A_n(1);
            obj.AtA(2,2) = obj.AtA(2,2) + A_n(2)*A_n(2);
            
            % A'*y += y_n*A_n';
            obj.AtY(1) = obj.AtY(1) + y_n*A_n(1);
            obj.AtY(2) = obj.AtY(2) + y_n*A_n(2);
        end       
        
        function x_est = ComputeEstimate(obj)
            % Compute the inverse of A'*A
            det = 1 / (obj.AtA(1,1)*obj.AtA(2,2) - obj.AtA(1,2)*obj.AtA(2,1));
            obj.AtA_inv = zeros(2,2);
            obj.AtA_inv(1,1) = det * obj.AtA(2,2);
            obj.AtA_inv(1,2) = -det * obj.AtA(1,2);
            obj.AtA_inv(2,1) = -det * obj.AtA(2,1);
            obj.AtA_inv(2,2) = det * obj.AtA(1,1);            
            
            % Compute estimate, AtA_inv * obj.AtY
            x_est = zeros(size(obj.AtA, 1), 1);
            x_est(1) = obj.AtA_inv(1,1)*obj.AtY(1) + obj.AtA_inv(1,2)*obj.AtY(2);
            x_est(2) = obj.AtA_inv(2,1)*obj.AtY(1) + obj.AtA_inv(2,2)*obj.AtY(2);            
        end
    end
end
