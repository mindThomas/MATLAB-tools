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
            z_hat = h + H * obj.x;
            % Compute innovation
            innov = z - z_hat;
            % Compute innovation variance
            S = H * obj.P * H' + R;
            % Compute Kalman gain
            K = obj.P * H' / S;
            % Compute posterior estimate
            obj.x = obj.x + K * innov;
            % Compute posterior estimation error variance
            %obj.P = obj.P - K * S * K';
            obj.P = obj.P - K * H * obj.P;
        end        
       
        function obj = dynamics_update(obj, y_new, f, g, Fx, Gx, Cov_qx, Cov_qy, Cov_qxy)
            % See http://user.it.uu.se/~thosc112/pubpdf/schongk2011.pdf
            % equation 11a to 12c
            % Discrete prediction model on the form:
            %   x[k] = f(y[k-1]) + Fx(y[k-1])*x[k-1] + Mx*qx[k-1]
            %   where qy[k] ~ N(0, Qy)                
            % Coupled with a non-linear state prediction model:
            %   y[k] = g(y[k-1]) + Gx(y[k-1]) * x[k-1] + My*qy[k-1]
            %   where qy[k] ~ N(0, Qy)                
            % With joint covariance (including correlation) given by:
            % q = [qx; qy]
            % Cov(q) = [Cov_qx  , Cov_qxy
            %           Cov_qxy', Cov_qy ]
            % This function handles how to incorporate a change in the
            % non-linear state into the linear state
            %    z[k-1] = y[k] - g(y[k-1])
            z = y_new - g;            
            
            % We assume additive Gaussian noise without any pre-mapping
            Mx = eye(size(Cov_qx));
            My = eye(size(Cov_qy));
            
            % Construct dynamics update transition matrix which takes the
            % correlation between qx and qy into account
            A = Fx - Mx * Cov_qxy' * inv(My*Cov_qy) * Gx;
            Q = Cov_qx  - Cov_qxy' * inv(Cov_qy) * Cov_qxy;
            
            % Compute innovation variance
            S = Gx * obj.P * Gx' + My * Cov_qy * My';
            
            % Compute Kalman gain
            K = Fx * obj.P * Gx' / S;
            
            % Compute posterior estimate
            obj.x = f + A * obj.x ...
                  + Mx * Cov_qxy' * inv(My*Cov_qy) * z ...
                  + K * (z - Gx*obj.x);
            
            % Compute posterior estimation error variance
            obj.P = Fx * obj.P * Fx' + Mx * Cov_qx * Mx' ...
                  - K * S * K';
        end             
        
    end
end