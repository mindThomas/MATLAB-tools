classdef RecursiveLeastSquares_KF
    properties
        X
        P
    end
    properties (Access = private)  
        Q
        R
        A_noise
    end    
    
    methods
        function obj = RecursiveLeastSquares_KF(X_init, P_init, Q, R, A_noise)            
            obj.X = X_init;
            obj.P = P_init;
            obj.Q = Q;
            obj.R = R;
            obj.A_noise = A_noise;
        end   
        
        function obj = Predict(obj)            
            %X_apriori = obj.X;
    
            % Determine model Jacobian (F)   -  OBS. This is not supposed to use/depend on apriori states!
            % Constant state model
            F_prev = eye(length(obj.X));    
    
            % Calculate apriori covariance of estimate error
            obj.P = F_prev * obj.P * F_prev' + obj.Q;     
        end
        
        function obj = Update(obj, A_n, y_n)
            % y_n = A_n*X
            z = y_n;
            
            % Compute measurement estimate
            z_hat = A_n * obj.X;            
            
            % Compute Jacobian
            H = A_n;
            
            % Calculate Kalman gain
            S = H * obj.P * H' + obj.X' * obj.A_noise * obj.X + obj.R;        
            K = obj.P * H' / S;  % P_apriori * H' * inv(S)        
 
            % Compute aposteriori
            obj.X = obj.X + K * (z - z_hat);    
            obj.P = (eye(length(obj.X)) - K*H) * obj.P;  
        end
    end
end
