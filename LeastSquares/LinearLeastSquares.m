classdef LinearLeastSquares
    % Different implementations of ordinary least squares
    % min ||A*x - y||^2
    % where || is the 2-norm (Euclidean distance)
    methods (Static)        
        function x = MATLAB(A, y)
            % MATLAB lscov based
            x = lscov(A, y);
        end
      
        function x = PseudoInverse(A, y)
            % Pseudo inverse based
            Ainv = inv(A'*A) * A';
            x = Ainv * y;
        end        
        
        function x = QR(A, y)
            % QR based
            % min ||A*x - y||^2
            [Q,R] = qr(A);
            QA = Q' * A;
            Qy = Q' * y;
            R_prime = QA(1:size(R,2),:);
            d = Qy(1:size(R,2),:);
            e = Qy(size(R,2):end,:);
            x = inv(R_prime) * d;
        end     
        
        function x = Cholesky(A, y)
            % Cholesky based
            % min ||A*x - y||^2
            R = chol(A'*A, 'upper');
            z = inv(R') * A'*y;
            x = inv(R) * z;
        end 
        
        function x = SVD(A, y)
            % SVD based
            % min ||A*x - y||^2
            [U,S,V] = svd(A);            
            singular_values_inv = 1./diag(S);
            singular_values_inv(find(min(diag(S)) <= eps)) = 0; % force the inverse of zero singular values to zero            
            S_inv = diag(singular_values_inv);
            % Since there are multiple solutions to the equation A*x = y we
            % chose the one with the smallest Euclidean norm, min |x|^2
            z = (U'*y);
            x = V * S_inv * z(1:length(singular_values_inv));
        end
        
        function x = SVD_homogeneous(A)
            % Solving an overdetermined system of Homogeneous Equations
            % min ||A*x||^2 
            % SVD based    
            % The solution x is the eigenvector corresponding to the only zero eigenvalue of A'*A
            % Equivalently, we can take the SVD of A and x is the column of V corresponding to the zero singular value of A
            [U,S,V] = svd(A);
            % Since the columns are ordered, this is the rightmost column of V
            x = V(:,end);
        end            
        
        function x = EIG_homogeneous(A)
            % Solving an overdetermined system of Homogeneous Equations
            % min ||A*x||^2 
            % Eigenvector decomposition based    
            % The solution x is the eigenvector corresponding to the only zero eigenvalue of A'*A            
            [V,D] = eig(A'*A);
            % Since the columns are ordered, this is the leftmost column of V
            x = V(:,1);
        end                      
    end
end
