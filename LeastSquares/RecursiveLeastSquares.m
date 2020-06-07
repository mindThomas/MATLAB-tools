classdef RecursiveLeastSquares
    % Different implementations of ordinary recursive least squares
    % min ||A*x - y||^2
    % where || is the 2-norm (Euclidean distance)
    % with an recursively incoming row of A_n and measurement y_n
    properties
        x_est
    end
    properties (Access = private)
        AtA
        AtA_inv
        AtY
    end    
    
    methods
        function obj = RecursiveLeastSquares(A_init, y_init)            
            % Pseudo inverse is given as
            % Ainv = inv(A'*A) * A';
            % x = Ainv * y;
            obj.AtA = A_init' * A_init;
            obj.AtY = A_init' * y_init;
            obj.AtA_inv = inv(obj.AtA);
            obj.x_est = zeros(size(A_init,2), 1);
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
            obj.AtA = obj.AtA + A_n'*A_n;
            obj.AtY = obj.AtY + y_n*A_n';
        end
        
        function obj = AddMeasurementAndComputeEstimate(obj, A_n, y_n)
            % A_n should be the new row of A to append
            if (size(A_n,1) ~= 1)
                error('A_n should be the new row of A to append');
            end
            if (size(A_n,2) ~= length(obj.AtY))
                error('Size of A_n does not match number of parameters to estimate');
            end
            % http://pfister.ee.duke.edu/courses/ece586/ex_proj_2008.pdf            
            obj.AtA = obj.AtA + A_n'*A_n;
            AtA_inv = inv(obj.AtA);
            
            % But this step can be computationally optimized by using the Sherman-Morrison formula to append the row and AtA
            AtA_inv2 = obj.AtA_inv - 1/(1 + A_n*obj.AtA_inv*A_n') * obj.AtA_inv*A_n'*A_n*obj.AtA_inv;            
            obj.AtA_inv = AtA_inv2; % store the updated inv(A'*A)
            
            % Use normal linear algebra to extend A'*y
            obj.AtY = obj.AtY + y_n*A_n';
            
            % Compute estimate
            obj.x_est = obj.AtA_inv * obj.AtY;
        end
        
        function x_est = ComputeEstimate(obj)
            % Compute the inverse of A'*A
            AtA_inv = inv(obj.AtA);
            
            % Compute estimate
            x_est = AtA_inv * obj.AtY;            
        end
    end
end
