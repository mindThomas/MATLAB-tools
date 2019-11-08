classdef KF
	properties %(SetAccess = private)
        x        
    end    
	properties (SetAccess = private)        
        P
    end        
    properties (SetAccess = private)%(Access = private)
        A
        B
        H
        Q
        R        
        x0
        u0   
        h0
        f0   
        Bq
        Hr
        Q_sqrt
        R_sqrt       
    end
    methods
        function obj = KF(varargin)                        

        end
        
        % Most advanced and complete function involving internal
        % differentiation (Jacobian) and discretization
        function obj = init_continuous_function(obj, ...
                                f, x0, u0, q0, ts, Qc, ... % process model
                                h, R, r0, ...              % measurement model
                                P0)                    % initial covariance
            % f is functions defined with f = @(x, u, q) ...            
            % h is functions defined with h = @(x, r) ...            
            if (max(abs(q0)) > eps)
                error('Process Noise is assumed to be zero-mean but input was different');
            end
            if (max(abs(r0)) > eps)
                error('Measurement Noise is assumed to be zero-mean but input was different');
            end         
            if (size(P0,1) ~= size(x0,1) || size(P0,2) ~= size(x0,1))
                error('Incorrect initial estimation error covariance matrix size');
            end
            
            % Prediction model is on the form:
            %   dx = f(x, u, q)          [ defined with dx = @(x, u, q) ... ]
            %   where q ~ N(0, Qc)
            % Which is first linearized into
            %   dx = f(x0, u0, 0) + Ac*(x-x0) + Bc*(u-u0) + Bq*q
            % where we have:
            %   Ac = df(x,u,q)/dx|x=x0,u=u0,q=0
            %   Bc = df(x,u,q)/du|x=x0,u=u0,q=0            
            %   Bq = df(x,u,q)/dq|x=x0,u=u0,q=0            
            [dfdx, dfdu, dfdq] = numjacobian3(f, x0, u0, q0);
            Ac = dfdx;
            Bc = dfdu;
            Bq = dfdq;
            % We compute the continuous noise covariance affecting dx as if
            %   dx = f(x, u, q) = f(x, u, 0) + Bq*q
            %   where q ~ N(0, Qc)
            Qcx = Bq * Qc * Bq';
            
            % Thereafter this is discretized into
            %   x[k+1] = x0 +
            %            ts * f(x0, u0) +
            %            Ad * (x[k]-x0) +
            %            Bd * (u[k]-u0) +
            %            q[k]
            %   where q[k] ~ N(0, ts*Qcx)    [ see "5.4.1 Selecting the 
            %   discrete time motion noise covariance" in ChM015x Course ]
            % Using the Euler method we have
            %   Ad = I + T*Ac
            %   Bd = T*Bc
            % Whereof the analytic method is:
            %   Ad = exp(Ac*ts)     [ = expm(Ac*ts) ]
            %   Bd = int_0^ts exp(A*tau) dtau * Bc
            % The integral can be computed from the Taylor expansion of exp(A*t)
            %   exp(A*t) = I + A*t + A^2*t^2/2 + A^3*t^3/3! + A^4*t^4/4! + ...
            % such that
            %   int_0^ts exp(A*tau) dtau = I*ts + A*ts^2/2 + A^2*ts^3/3! + ...
            [Ad, Bd] = discretize_analytic(Ac, Bc, ts);       
            
            % Measurement model is on the form:
            %   z = h(x, r)          [ defined with h = @(x,r) ... ]
            %   where r ~ N(0, R)
            % Which is linearized into
            %   z = h(x0, 0) + Hx*(x-x0) + Hr*r
            % where:
            %   Hx = dh(x,r)/dx|x=x0,r=0
            %   Hr = dh(x,r)/dr|x=x0,r=0
            [dhdx, dhdr] = numjacobian2(h, x0, r0);  
            % We compute the noise covariance affecting dx as if
            %   z = h(x, 0) + Hr*r
            %   where r ~ N(0, R)
            Rx = dhdr * R * dhdr';
 
            obj.A = Ad;
            obj.B = Bd;
            obj.H = dhdx;
            obj.Q = ts * Qcx;
            obj.R = Rx;
            obj.x0 = x0;
            obj.u0 = u0;
            obj.h0 = h(x0,r0);
            obj.f0 = ts * f(x0,u0,q0);        
            obj.Bq = Bq;
            obj.Hr = dhdr;
            obj.Q_sqrt = chol(Q, 'lower');
            obj.R_sqrt = chol(R, 'lower');
            obj.x = x0;
            obj.P = P0;
        end          
        
        % Initialize Kalman filter object based on process and measurement
        % model functions and Jacobian functions        
        function obj = init_continuous_function_jacobians(obj, ...
                            f, Fx, Fu, Fq, x0, u0, q0, ts, Qc, ... % process model
                            h, Hx, Hr, R, r0, ...                   % measurement model
                            P0)                            % initial covariance
            % f, Fx, Fu, Fq is functions defined with f = @(x, u, q) ...
            % h, H is functions defined with h = @(x, r) ...
            if (max(abs(q0)) > eps)
                error('Process Noise is assumed to be zero-mean but input was different');
            end
            if (max(abs(r0)) > eps)
                error('Measurement Noise is assumed to be zero-mean but input was different');
            end    
            if (size(P0,1) ~= size(x0,1) || size(P0,2) ~= size(x0,1))
                error('Incorrect initial estimation error covariance matrix size');
            end
            
            % Prediction model is on the form:
            %   dx = f(x, u, q)
            %   where q ~ N(0, Qc)
            % Which is first linearized into
            %   dx = f(x0, u0, 0) + Ac*(x-x0) + Bc*(u-u0) + Bq*q
            Ac = Fx(x0,u0,q0);
            Bc = Fu(x0,u0,q0);
            Bq = Fq(x0,u0,q0);
            % We compute the continuous noise covariance affecting dx as if
            %   dx = f(x, u, q) = f(x, u, 0) + Bq*q
            %   where q ~ N(0, Qc)
            Qcx = Bq * Qc * Bq';    
            
            % Thereafter this is discretized into
            %   x[k+1] = x0 +
            %            ts * f(x0, u0) +
            %            Ad * (x[k]-x0) +
            %            Bd * (u[k]-u0) +
            %            q[k]
            %   where q[k] ~ N(0, ts*Qcx)    [ see "5.4.1 Selecting the 
            %   discrete time motion noise covariance" in ChM015x Course ]
            % Using the Euler method we have
            %   Ad = I + T*Ac
            %   Bd = T*Bc
            % Whereof the analytic method is:
            %   Ad = exp(Ac*ts)     [ = expm(Ac*ts) ]
            %   Bd = int_0^ts exp(A*tau) dtau * Bc
            % The integral can be computed from the Taylor expansion of exp(A*t)
            %   exp(A*t) = I + A*t + A^2*t^2/2 + A^3*t^3/3! + A^4*t^4/4! + ...
            % such that
            %   int_0^ts exp(A*tau) dtau = I*ts + A*ts^2/2 + A^2*ts^3/3! + ...
            [Ad, Bd] = discretize_analytic(Ac, Bc, ts);    
            
            % Measurement model is on the form:
            %   z = h(x, r)          [ defined with h = @(x,r) ... ]
            %   where r ~ N(0, R)
            % Which is linearized into
            %   z = h(x0, 0) + Hx*(x-x0) + Hr*r
            % where:
            %   Hx = dh(x,r)/dx|x=x0,r=0
            %   Hr = dh(x,r)/dr|x=x0,r=0              
            % We compute the noise covariance affecting dx as if
            %   z = h(x, 0) + Hr*r
            %   where r ~ N(0, R)
            Rx = Hr(x0,r0) * R * Hr(x0,r0)';            
 
            obj.A = Ad;
            obj.B = Bd;
            obj.H = Hx(x0,r0);
            obj.Q = ts * Qcx;
            obj.R = Rx;
            obj.x0 = x0;
            obj.u0 = u0;
            obj.h0 = h(x0,r0);
            obj.f0 = ts * f(x0,u0,q0);   
            obj.Bq = Bq;
            obj.Hr = Hr(x0,r0);
            obj.Q_sqrt = chol(Q, 'lower');
            obj.R_sqrt = chol(R, 'lower');
            obj.x = x0;
            obj.P = P0;
        end                
                
        function obj = init_discrete_function(obj, ...
                                f, x0, u0, q0, Q, ...      % process model
                                h, R, r0, ...              % measurement model
                                P0)                    % initial covariance
            % f is defined as "x[k+1] = @(x[k], u[k], q[k]) ... "
            % h is defined as "z[k] = @(x[k], r[k]) ... "
            if (max(abs(q0)) > eps)
                error('Process Noise is assumed to be zero-mean but input was different');
            end
            if (max(abs(r0)) > eps)
                error('Measurement Noise is assumed to be zero-mean but input was different');
            end    
            if (size(P0,1) ~= size(x0,1) || size(P0,2) ~= size(x0,1))
                error('Incorrect initial estimation error covariance matrix size');
            end            
            
            % Prediction model is on the form:
            %   x[k+1] = f(x[k], u[k], q[k])
            % Where q[k] ~ N(0, Q)            
            
            % Which is linearized into
            %   x[k+1] = f(x0, u0, 0) + Ad*(x-x0) + Bd*(u-u0) + Bq*q
            % where we have:
            %   Ad = df(x,u)/dx|x=x0,u=u0,q=0
            %   Bd = df(x,u)/du|x=x0,u=u0,q=0            
            %   Bq = df(x,u)/du|x=x0,u=u0,q=0            
            [dfdx, dfdu, dfdq] = numjacobian3(f, x0, u0, q0);
            Ad = dfdx;
            Bd = dfdu;
            Bq = dfdq;
            % We compute the discrete noise covariance affecting x[k] as if
            %   x[k] = f(x[k], u[k], q[k]) = f(x[k], u[k], 0) + Bq*q[k]
            %   where q ~ N(0, Qc)
            Qx = Bq * Q * Bq';    
            
            % Measurement model is on the form:
            %   z[k] = h(x[k], r[k])
            %   where r[k] ~ N(0, R)
            % Which is linearized into
            %   z = h(x0) + H*(x-x0) + Hr*r
            % where:
            %   Hx = dh(x)/dx|x=x0,r=0         
            %   Hr = dh(x)/dr|x=x0,r=0         
            [dhdx, dhdr] = numjacobian2(h, x0, r0);
            % We compute the noise covariance affecting dx as if
            %   z = h(x, 0) + Hr*r
            %   where r ~ N(0, R)
            Rx = dhdr * R * dhdr'; 
 
            obj.A = Ad;
            obj.B = Bd;
            obj.H = dhdx;
            obj.Q = Qx;
            obj.R = Rx;
            obj.x0 = x0;
            obj.u0 = u0;
            obj.h0 = h(x0);
            obj.f0 = f(x0,u0,q0) - x0; 
            obj.Bq = Bq;
            obj.Hr = dhdr(x0,r0);            
            obj.Q_sqrt = chol(Q, 'lower');
            obj.R_sqrt = chol(R, 'lower');
            obj.x = x0;
            obj.P = P0;                       
        end
        
        function obj = init_discrete_function_jacobians(obj, ...
                                f, Fx, Fu, Fq, x0, u0, q0, Q, ... % process model
                                h, Hx, Hr, R, r0, ...              % measurement model
                                P0)                       % initial covariance
            % f is defined as "x[k+1] = @(x[k], u[k], q[k]) ... "
            % h is defined as "z[k] = @(x[k], r[k]) ... "
            if (max(abs(q0)) > eps)
                error('Process Noise is assumed to be zero-mean but input was different');
            end
            if (max(abs(r0)) > eps)
                error('Measurement Noise is assumed to be zero-mean but input was different');
            end   
            if (size(P0,1) ~= size(x0,1) || size(P0,2) ~= size(x0,1))
                error('Incorrect initial estimation error covariance matrix size');
            end            
            
            % Prediction model is on the form:
            %   x[k+1] = f(x[k], u[k], q[k])
            % Where q[k] ~ N(0, Q)            
            
            % Which is linearized into
            %   x[k+1] = f(x0, u0, 0) + Ad*(x-x0) + Bd*(u-u0) + Bq*q
            Ad = Fx(x0, u0, q0);
            Bd = Fu(x0, u0, q0);
            Bq = Fq(x0, u0, q0);
            % We compute the discrete noise covariance affecting x[k] as if
            %   x[k] = f(x[k], u[k], q[k]) = f(x[k], u[k], 0) + Bq*q[k]
            %   where q ~ N(0, Qc)
            Qx = Bq * Q * Bq';   

            % Measurement model is on the form:
            %   z[k] = h(x[k], r[k])
            %   where r[k] ~ N(0, R)
            % Which is linearized into
            %   z = h(x0) + H*(x-x0) + Hr*r
            % where:
            %   Hx = dh(x)/dx|x=x0,r=0         
            %   Hr = dh(x)/dr|x=x0,r=0         
            % We compute the noise covariance affecting dx as if
            %   z = h(x, 0) + Hr*r
            %   where r ~ N(0, R)
            Rx = Hr(x0,q0) * R * Hr(x0,q0)';                     
 
            obj.A = Ad;
            obj.B = Bd;
            obj.H = Hx(x0,r0);
            obj.Q = Qx;
            obj.R = Rx;
            obj.x0 = x0;
            obj.u0 = u0;
            obj.h0 = h(x0,r0);
            obj.f0 = f(x0,u0,q0) - x0; 
            obj.Bq = Bq;
            obj.Hr = Hr(x0,r0);               
            obj.Q_sqrt = chol(Q, 'lower');
            obj.R_sqrt = chol(R, 'lower');
            obj.x = x0;
            obj.P = P0;                       
        end        
        
        % Simplest form of the Kalman filter with x0=0 and u0=0
        function obj = init_discrete(obj, A, B, H, Q, R, P0)
            % Prediction model is on the form:
            %   x[k+1] = A*x[k] + B*u[k] + q[k]
            %   where q[k] ~ N(0, Q)
            % Measurement model is on the form:
            %   z[k] = H*x[k] + r[k]
            %   where r[k] ~ N(0, R)
            obj.A = A;
            obj.B = B;
            obj.H = H;
            obj.Q = Q;
            obj.R = R;
            obj.x0 = zeros(size(A,1), 1);
            obj.u0 = zeros(size(B,1), 1);
            obj.h0 = zeros(size(H,1), 1);
            obj.f0 = f(x0,u0,q0) - x0;   
            obj.Bq = eye(size(Q));
            obj.Hr = eye(size(R));               
            obj.Q_sqrt = chol(Q, 'lower');
            obj.R_sqrt = chol(R, 'lower');
            obj.x = obj.x0;
            obj.P = P0;
        end

                        
        function obj = predict(obj, u) 
            % Compute prior estimate
            obj.x = obj.x0 + ...
                    obj.f0 + ...
                    obj.A * (obj.x - obj.x0);
            if (~isempty(obj.B))
                obj.x = obj.x + ...
                        obj.B * (u - obj.u0);                    
            end
            % Compute prior estimation error variance
            obj.P = obj.A * obj.P * obj.A' + obj.Q;
        end
        
        function obj = update(obj, z)
            z_hat = obj.h0 + obj.H * (obj.x - obj.x0);
            % Compute innovation
            innov = z - z_hat;
            % Compute innovation variance
            S = obj.H * obj.P * obj.H' + obj.R;
            % Compute Kalman gain
            K = obj.P * obj.H' * inv(S);
            % Compute posterior estimate
            obj.x = obj.x + K * innov;
            % Compute posterior estimation error variance
            obj.P = obj.P - K * S * K';
        end
    
        function x = sampleProcessModel(obj, u)
            % Generate a realization from the process model based on x and Q
            % Generate noise realization, q[k] ~ N(0, Q)
            q = obj.Q_sqrt * randn([size(obj.Bq), 2]);
            % Apply the process model
            x = obj.x0 + ...
                obj.f0 + ...
                obj.A * (obj.x - obj.x0) + ...                
                obj.Bq * q;   
            if (~isempty(obj.B))
                x = x + ...
                    obj.B * (u - obj.u0);
            end            
        end 
        
        function z = sampleMeasurementModel(obj)
            % Generate a realization from the measurement model based on x and R
            % Generate noise realization, r ~ N(0, R)
            r = obj.R_sqrt * randn([size(obj.Hr), 2]);
            % Apply the linearized measurement model
            z = obj.h0 + obj.H * (obj.x - obj.x0) + obj.Hr * r;
        end 
        
        function x = samplePosterior(obj)
           % Generate a realization of the state based on x and P
           P_sqrt = chol(obj.P, 'lower');
           % Generate estimation error noise realization
           p = P_sqrt * randn([length(obj.P), 1]);
           % Compose estimate and noise into realization
           x = obj.x + p;
        end
                
        function z = sampleMeasurement(obj)
            % Generate a realization of a measurement based on a
            % realization of the state and the measurement noise
            x = samplePosterior();
            % Generate noise realization, r ~ N(0, R)
            r = obj.R_sqrt * randn([size(obj.Hr), 2]);
            % Apply the linearized measurement model
            z = obj.h0 + obj.H * (x - obj.x0) + obj.Hr * r;
        end        
    end
end