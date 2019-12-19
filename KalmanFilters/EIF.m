classdef EIF
    % Extended Information Filter
    % The information filter implemented the Kalman filter in Logarithmic form
    % The key difference between the EKF and the EIF arises from the way
    % the Gaussian belief is represented. Whereas in the Kalman filter family of algorithms,
    % Gaussians are represented by their moments (mean, covariance), information filters
    % represent Gaussians in their canonical representation, which is comprised of an information
    % matrix and an information vector.
    % 
    % See Chapter 3.4.3 of Probabilistic Robotics by Sebastian Thrun
    % See also http://ais.informatik.uni-freiburg.de/teaching/ws12/mapping/pdf/slam06-eif.pdf
    %
    % The measurement update is the difficult step in Kalman filters.
    % But for the information filter measurement updates are additive and
    % are even more efficient if measurements carry only information about a subset of
    % all state variables at a time.
    %
    % Note how the information filter supports global/infinite uncertainty
    % by setting the information matrix to 0
    %
    % OBS: For high dimensional state spaces, the information filter is
    % generally believed to be computationally inferior to (worse than) the Kalman filter.
    % In fact, this is one of the reasons why the EKF has been vastly more
    % popular than the EIF.
	properties %(SetAccess = private)
        xi % information vector = P^-1 * mu = Omega * mu
    end    
	properties (SetAccess = private)        
        Omega % information matrix = is the inverse of the covariance matrix P
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
        function obj = EIF(varargin)                        

        end        
        
        function [x, P] = getEstimate(obj)
            % Function to convert from information vector and
            % information matrix to estimate mean and covarinace
            P = inv(obj.Omega);
            x = P * obj.xi;
        end

        function obj = setState(obj, x)
            % Function to convert from estimate mean to information vector
            obj.xi = obj.Omega * x;
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
            obj.Omega = inv(P_init);
            obj.xi = obj.Omega * x_init;            
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
            obj.Omega = inv(P_init);
            obj.xi = obj.Omega * x_init;                        
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
            obj.Omega = inv(P_init);
            obj.xi = obj.Omega * x_init;            
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
            obj.Omega = inv(P_init);
            obj.xi = obj.Omega * x_init;            
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
            % Compute current state
            P = inv(obj.Omega);
            x = P * obj.xi;
            
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
                dhdx = obj.Hx(x, r);
                dhdr = obj.Hr(x, r);
            else
                [dhdx, dhdr] = numjacobian2(obj.h, x, r);  
            end            
            H = dhdx;
            % The above is written differently with noise affecting the states directly
            %   z = h(x_prior, 0) + Hx*(x-x_prior) + r_x
            %   where r_x ~ N(0, Hr*R*Hr')            
            Rx = dhdr * obj.R * dhdr';                                   
            
            % Compute posterior information matrix
            invR = inv(Rx);
            obj.Omega = obj.Omega + H' * invR * H;
            
            % Compute predicted measurement
            z_hat = obj.h(x, r); % with r=0
            % Compute posterior information vector
            obj.xi = obj.xi + H' * invR * (z - z_hat + H * x);    
            % Note that this equation is different from Table 3.5 in
            % Probabilistic Robotics, which has an sign in front of H*x. Instead
            % the correct equation is shown in http://ais.informatik.uni-freiburg.de/teaching/ws12/mapping/pdf/slam06-eif.pdf
        end             
        
    end
    methods (Access = private)
        function obj = predict_discrete(obj, varargin)
            if (nargin == 2)
                u = varargin{1};
            else
                u = [];
            end
            % Compute previous state
            P = inv(obj.Omega);
            x_prev = P * obj.xi;
            
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
                Ad = obj.Fx(x_prev,u,q);
                %Bd = obj.Fu(x_prev,u,q);
                Bdq = obj.Fq(x_prev,u,q);
            else
                % Perform numerical differentiation to compute Jacobians numerically
                [dfdx, dfdu, dfdq] = numjacobian3(obj.f, x_prev, u, q);
                Ad = dfdx;
                %Bd = dfdu;
                Bdq = dfdq;
            end 
            
            % Compute prior information matrix            
            P_prior = Ad * P * Ad' + Bdq * obj.Q * Bdq';
            obj.Omega = inv(P_prior);            
                        
            % Compute prior estimate using the discrete prediction model:
            %   x[k+1] = f(x_hat[k], u[k], 0)            
            x = obj.f(x_prev, u, q);
            obj.xi = obj.Omega * x;           
        end       
        
        function obj = predict_continuous(obj, varargin)
            if (nargin == 2)
                u = varargin{1};
            else
                u = [];
            end
            % Compute previous state
            P = inv(obj.Omega);
            x_prev = P * obj.xi;
            
            % Continuous prediction model on the form:
            %   dx = f(x, u, q)          [ defined with dx = @(x, u, q) ... ]
            %   where q ~ N(0, Qc)
            % We linearize the system with respect to noise such that
            %   dx = f(x, u, q) = f(x, u, 0) + Bq*q
            q = zeros(size(obj.Q,1),1); % no noise
            if (obj.JacobiansAvailable)
                % Evaluate Jacobians at previous posterior state
                Ac = obj.Fx(x_prev,u,q);
                Bc = obj.Fu(x_prev,u,q);
                Bq = obj.Fq(x_prev,u,q);
            else
                % Perform numerical differentiation to compute Jacobians numerically
                [dfdx, dfdu, dfdq] = numjacobian3(obj.f, x_prev, u, q);
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
            
            % Compute prior information matrix            
            P_prior = Ad * P * Ad' + Qd;
            obj.Omega = inv(P_prior);            
                        
            % Compute prior estimate using the continuous prediction model:
            %   dx = f(x, u, 0)
            % OBS: This can be done in 3 ways:
            %   1. Proper integration from t=t0 to t=t0+ts (e.g. with ODE solver)
            %   2. Euler method as listed above            
            % Forward Euler method
            %   x[k+1] = x[k] + ts*f(x[k], u[k], 0)            
            x = x_prev + obj.ts*obj.f(x_prev, u, q); % with q=0
            obj.xi = obj.Omega * x;              
        end                    
    end
end