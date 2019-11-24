classdef RaoBlackwellizedParticle
	properties %(SetAccess = private)        
        weight
        nonlinear_state
        prev_nonlinear_state
        % The linear state(s) is/are represented with a Gaussian and are
        % filtered using a linear Kalman filter
        x
        P
    end        
    methods
        function obj = RaoBlackwellizedParticle(weight, nonlinear_state, linear_state, linear_covariance)
            obj.weight = weight;
            obj.nonlinear_state = nonlinear_state;
            obj.prev_nonlinear_state = nonlinear_state;
            obj.x = linear_state;
            obj.P = linear_covariance;            
        end     
        
        function obj = linear_predict(obj, f, A, Q)
            % Discrete prediction model on the form:
            %   x[k+1] = f + A*x[k] + q[k]
            %   where q[k] ~ N(0, Q)                                               
            % Compute prior estimate using the discrete prediction model:
            %   x_hat[k+1] = A*x_hat[k]
            obj.x = f + A * obj.x;
            % Compute prior estimation error variance
            obj.P = A * obj.P * A' + Q; 
        end  
        
        function obj = linear_update(obj, h, H, R, z)
            % Measurement model is on the form:
            %   z = h + H*x + r
            %   where r ~ N(0, R)                              
            % From the measurement model a predicted measurement is computed
            z_hat = H * obj.x;
            % Compute innovation
            innov = z - z_hat;
            % Compute innovation variance
            S = H * obj.P * H' + R;
            % Compute Kalman gain
            K = obj.P * H' / S;
            % Compute posterior estimate
            obj.x = obj.x + K * innov;
            % Compute posterior estimation error variance
            %obj.P = obj.P - K * S * obj.K';
            obj.P = obj.P - K * H * obj.P;
        end        
        
    end
end