classdef EKF
	properties %(SetAccess = private)
        x        
    end    
	properties (SetAccess = private)        
        P
        K
        S
    end        
    properties (SetAccess = private)%(Access = private)
        f        
        h
        Q
        R           
        Q_sqrt
        R_sqrt 
        ts
    end
    properties (SetAccess = private)%(Access = private)
        JacobiansAvailable
        Fx
        Fu
        Fq
        Hx
        Hr
    end    
    methods
        function obj = EKF(varargin)                        

        end        

        function obj = init_continuous(obj, f, Qc, ts, h, R, x_init, P_init)
            % Prediction model is on the form:
            %   dx = f(x, u, q)
            %   where q ~ N(0, Q)
            % Measurement model is on the form:
            %   z = h(x, r)
            obj.f = f;
            obj.h = h;
            obj.JacobiansAvailable = false;
            obj.x = x_init;
            obj.P = P_init; 
            obj.ts = ts;
            obj.Q = Qc;
            obj.R = R;
            obj.Q_sqrt = chol(Qc, 'lower');
            obj.R_sqrt = chol(R, 'lower');            
        end
        
        function obj = init_continuous_jacobians(obj, f, Fx, Fu, Fq, Qc, ts, h, Hx, Hr, R, x_init, P_init)
            % Prediction model is on the form:
            %   dx = f(x, u, q)
            %   where q ~ N(0, Q)
            % Measurement model is on the form:
            %   z = h(x, r)
            obj.f = f;
            obj.h = h;
            obj.Fx = Fx;
            obj.Fu = Fu;
            obj.Fq = Fq;
            obj.Hx = Hx;
            obj.Hr = Hr;            
            obj.JacobiansAvailable = true;
            obj.x = x_init;
            obj.P = P_init; 
            obj.ts = ts;
            obj.Q = Qc;
            obj.R = R;
            obj.Q_sqrt = chol(Qc, 'lower');
            obj.R_sqrt = chol(R, 'lower');            
        end
        
        function obj = init_discrete(obj, f, Q, h, R, x_init, P_init)
            % Prediction model is on the form:
            %   x[k+1] = f(x[k], u[k], q[k])
            %   where q[k] ~ N(0, Q)
            % Measurement model is on the form:
            %   z[k] = h(x[k], r[k])
            %   where r[k] ~ N(0, R)   
            obj.f = f;
            obj.h = h;      
            obj.JacobiansAvailable = false;
            obj.x = x_init;
            obj.P = P_init;
            obj.ts = 0;
            obj.Q = Q;
            obj.R = R;
            obj.Q_sqrt = chol(Q, 'lower');
            obj.R_sqrt = chol(R, 'lower');
        end
        
        function obj = init_discrete_jacobians(obj, f, Fx, Fu, Fq, Q, h, Hx, Hr, R, x_init, P_init)
            % Prediction model is on the form:
            %   x[k+1] = f(x[k], u[k], q[k])
            %   where q[k] ~ N(0, Q)
            % Measurement model is on the form:
            %   z[k] = h(x[k], r[k])
            %   where r[k] ~ N(0, R)   
            obj.f = f;
            obj.h = h;      
            obj.Fx = Fx;
            obj.Fu = Fu;
            obj.Fq = Fq;
            obj.Hx = Hx;
            obj.Hr = Hr;
            obj.JacobiansAvailable = true;
            obj.x = x_init;
            obj.P = P_init;
            obj.ts = 0;
            obj.Q = Q;
            obj.R = R;
            obj.Q_sqrt = chol(Q, 'lower');
            obj.R_sqrt = chol(R, 'lower');
        end
                        
        function obj = predict(obj, varargin) 
            if (nargin == 2)
                u = varargin{1};
            else
                u = [];
            end            
            if (obj.ts == 0)
                obj = predict_discrete(obj, u);                               
            else
                obj = predict_continuous(obj, u);
            end           
        end
        
        function obj = update(obj, z)
            % Measurement model is on the form:
            %   z = h(x, r)          [ defined with h = @(x,r) ... ]
            %   where r ~ N(0, R)
            % Which is linearized into
            %   z = h(x_prior, 0) + Hx*(x-x_prior) + Hr*r
            % where:
            %   Hx = dh(x,r)/dx|x=x0,r=0
            %   Hr = dh(x,r)/dr|x=x0,r=0
            r = zeros(size(obj.R,1),1); % no noise
            if (obj.JacobiansAvailable)
                dhdx = obj.Hx(obj.x, r);
                dhdr = obj.Hr(obj.x, r);
            else
                [dhdx, dhdr] = numjacobian2(obj.h, obj.x, r);  
            end            
            H = dhdx;
            % The above is written differently with noise affecting the states directly
            %   z = h(x_prior, 0) + Hx*(x-x_prior) + r_x
            %   where r_x ~ N(0, Hr*R*Hr')            
            Rx = dhdr * obj.R * dhdr';           
            
            % From the measurement model a predicted measurement is computed
            z_hat = obj.h(obj.x, r); % with r=0
            % Compute innovation
            innov = z - z_hat;
            % Compute innovation variance
            obj.S = H * obj.P * H' + Rx;
            % Compute Kalman gain
            obj.K = obj.P * H' / obj.S;
            % Compute posterior estimate
            obj.x = obj.x + obj.K * innov;
            % Compute posterior estimation error variance
            %obj.P = obj.P - obj.K * obj.S * obj.K';
            obj.P = obj.P - obj.K * H * obj.P;
        end
    
        function x = sampleProcessModel(obj, u)
            % Generate a realization from the process model based on x and Q
            % Generate noise realization, q[k] ~ N(0, Q)
            q = obj.Q_sqrt * randn([size(obj.Bq, 2), 1]);
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
            r = obj.R_sqrt * randn([size(obj.Hr, 2), 1]);
            % Apply the linearized measurement model
            z = obj.h0 + obj.H * (obj.x - obj.x0) + obj.Hr * r;
        end 
        
        function x = samplePosterior(obj)
           % Generate a realization of the state based on x and P
           P_sqrt = chol(obj.P, 'lower');
           % Generate estimation error noise realization
           p = P_sqrt * randn([size(obj.P, 1), 1]);
           % Compose estimate and noise into realization
           x = obj.x + p;
        end
                
        function z = sampleMeasurement(obj)
            % Generate a realization of a measurement based on a
            % realization of the state and the measurement noise
            x = samplePosterior();
            % Generate noise realization, r ~ N(0, R)
            r = obj.R_sqrt * randn([size(obj.Hr, 2), 1]);
            % Apply the linearized measurement model
            z = obj.h0 + obj.H * (x - obj.x0) + obj.Hr * r;
        end           
        
    end
    methods (Access = private)
        function obj = predict_discrete(obj, varargin)
            if (nargin == 2)
                u = varargin{1};
            else
                u = [];
            end
            % Discrete prediction model on the form:
            %   x[k+1] = f(x[k], u[k], q[k])
            %   where q[k] ~ N(0, Q)                        
            % Which is linearized into
            %   x[k+1] = f(x_hat[k], u[k], 0) + Ad*(x[k]-x_hat[k]) + Bd*u[k] + Bq*q[k]
            % where we have:
            %   Ad =  dx[k+1]/dx[k]|x=x_hat[k],u=u[k],q=0
            %   Bd =  dx[k+1]/du[k]|x=x_hat[k],u=u[k],q=0            
            %   Bdq = d x[k+1]/dq[k]|x=x_hat[k],u=u[k],q=0                 
            q = zeros(size(obj.Q,1),1); % no noise
            if (obj.JacobiansAvailable)
                % Evaluate Jacobians at previous posterior state
                Ad = obj.Fx(obj.x,u,q);
                %Bd = obj.Fu(obj.x,u,q);
                Bdq = obj.Fq(obj.x,u,q);
            else
                % Perform numerical differentiation to compute Jacobians numerically
                [dfdx, dfdu, dfdq] = numjacobian3(obj.f, obj.x, u, q);
                Ad = dfdx;
                %Bd = dfdu;
                Bdq = dfdq;
            end 
            
            % Compute prior estimate using the discrete prediction model:
            %   x[k+1] = f(x_hat[k], u[k], 0)
            obj.x = obj.f(obj.x, u, q);

            % Compute prior estimation error variance
            obj.P = Ad * obj.P * Ad' + Bdq * obj.Q * Bdq'; 
        end       
        
        function obj = predict_continuous(obj, varargin)
            if (nargin == 2)
                u = varargin{1};
            else
                u = [];
            end
            % Continuous prediction model on the form:
            %   dx = f(x, u, q)          [ defined with dx = @(x, u, q) ... ]
            %   where q ~ N(0, Qc)
            % We linearize the system with respect to noise such that
            %   dx = f(x, u, q) = f(x, u, 0) + Bq*q
            q = zeros(size(obj.Q,1),1); % no noise
            if (obj.JacobiansAvailable)
                % Evaluate Jacobians at previous posterior state
                Ac = obj.Fx(obj.x,u,q);
                Bc = obj.Fu(obj.x,u,q);
                Bq = obj.Fq(obj.x,u,q);
            else
                % Perform numerical differentiation to compute Jacobians numerically
                [dfdx, dfdu, dfdq] = numjacobian3(obj.f, obj.x, u, q);
                Ac = dfdx;
                Bc = dfdu;
                Bq = dfdq;
            end
            % This thus corresponds to 
            %   dx = f(x, u, 0) + q_x
            %   where q_x ~ N(0, Bq*Qc*Bq')            
            Qcx = Bq * obj.Q * Bq';           
                        
            % Next the continuous dynamics is discretized using Eulers method into
            %   x[k+1] = x[k] + ts*f(x[k], u[k], 0) + q[k]
            %   where q[k] ~ N(0, ts*Qcx)     [ see "5.4.1 Selecting the 
            %    discrete time motion noise covariance" in ChM015x Course ]   
            % Which is linearized into
            %   x[k+1] = x[k] + ts*f(x_hat[k], u[k], 0) + ts*Ac*(x[k]-x_hat[k]) + ts*Bc*u[k] + q[k]
            % where we have:
            %   Ac = df(x,u,q)/dx|x=x_hat,u=u,q=0
            %   Bc = df(x,u,q)/du|x=x_hat,u=u,q=0            
            %   Bq = df(x,u,q)/dq|x=x_hat,u=u,q=0 
            %Ad = obj.ts*Ac;
            %Bd = obj.ts*Bc;
            Qd = obj.ts*Qcx;            
            % OBS. The above could possible also have been replaced by
            % Matrix exponential discretization:
            [Ad, Bd] = discretize_analytic(Ac, Bc, obj.ts); 
            
            % Compute prior estimate using the continuous prediction model:
            %   dx = f(x, u, 0)
            % OBS: This can be done in 3 ways:
            %   1. Proper integration from t=t0 to t=t0+ts (e.g. with ODE solver)
            %   2. Euler method as listed above            
            % Forward Euler method
            %   x[k+1] = x[k] + ts*f(x[k], u[k], 0)
            obj.x = obj.x + ts*obj.f(obj.x, u, q); % with q=0

            % Compute prior estimation error variance
            obj.P = Ad * obj.P * Ad' + Qd;             
        end                    
    end
end