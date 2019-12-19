classdef IF
    % Information Filter
    % The information filter implemented the Kalman filter in Logarithmic form
    % The key difference between the KF and the IF arises from the way
    % the Gaussian belief is represented. Whereas in the Kalman filter family of algorithms,
    % Gaussians are represented by their moments (mean, covariance), information filters
    % represent Gaussians in their canonical representation, which is comprised of an information
    % matrix and an information vector.
    % 
    % See Chapter 3.4.3 of Probabilistic Robotics by Sebastian Thrun
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
	properties %(SetAccess = private)
        xi % information vector = P^-1 * mu = Omega * mu
    end    
	properties (SetAccess = private)        
        Omega % information matrix = is the inverse of the covariance matrix P
        K
        S
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
        function obj = IF(varargin)                        

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
            %Qcx = Bq * Qc * Bq';
            % NO. Instead of assuming the noise covariance to be ts*Qcx we
            % assume the noise to be applied and held throughout the sample
            % period, similar to an input. See below.           
            
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
            
            % We perform a similar discretization with regards to the noise input
            [Ad_, Bqd] = discretize_analytic(Ac, Bq, ts);  
            % The noise is then computed as
            Q = Bqd * Qc * Bqd';
            % This noise discretization can be improved to an analytically correct version as described in
            % http://webee.technion.ac.il/people/shimkin/Estimation09/ch8_target.pdf
            
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
            obj.Q = Q;
            obj.R = Rx;
            obj.x0 = x0;
            obj.u0 = u0;
            obj.h0 = h(x0,r0);
            obj.f0 = ts * f(x0,u0,q0);        
            obj.Bq = Bq;
            obj.Hr = dhdr;
            obj.Q_sqrt = chol(Qc, 'lower');
            obj.R_sqrt = chol(R, 'lower');            
            obj.Omega = inv(P0);
            obj.xi = obj.Omega * x0;
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
            %Qcx = Bq * Qc * Bq';
            % NO. Instead of assuming the noise covariance to be ts*Qcx we
            % assume the noise to be applied and held throughout the sample
            % period, similar to an input. See below.              
            
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
            
            % We perform a similar discretization with regards to the noise input
            [Ad_, Bqd] = discretize_analytic(Ac, Bq, ts);  
            % The noise is then computed as
            Q = Bqd * Qc * Bqd';
            % This noise discretization can be improved to an analytically correct version as described in
            % http://webee.technion.ac.il/people/shimkin/Estimation09/ch8_target.pdf
            
            % This noise discretization can be improved to an analytically correct version as described in
            % http://webee.technion.ac.il/people/shimkin/Estimation09/ch8_target.pdf
            
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
            obj.Q = Q;
            obj.R = Rx;
            obj.x0 = x0;
            obj.u0 = u0;
            obj.h0 = h(x0,r0);
            obj.f0 = ts * f(x0,u0,q0);   
            obj.Bq = Bq;
            obj.Hr = Hr(x0,r0);
            obj.Q_sqrt = chol(Qc, 'lower');
            obj.R_sqrt = chol(R, 'lower');
            obj.Omega = inv(P0);
            obj.xi = obj.Omega * x0;
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
            obj.Omega = inv(P0);
            obj.xi = obj.Omega * x0;
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
            %   where q[k] ~ N(0, Q)                        
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
            obj.Omega = inv(P0);
            obj.xi = obj.Omega * x0;
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
            obj.Omega = inv(P0);
            obj.xi = obj.Omega * obj.x0;
        end

                        
        function obj = predict(obj, varargin) 
            if (nargin == 2)
                u = varargin{1};
            else
                u = [];
            end
            % Compute prior information matrix
            P = inv(obj.Omega);
            P_prior = obj.A * P * obj.A' + obj.Q;
            obj.Omega = inv(P_prior);
            
            % Compute prior information vector
            x_prev = P * obj.xi;
            x_prior = obj.x0 + ...
                    obj.f0 + ...
                    obj.A * (x_prev - obj.x0);
            if (~isempty(obj.B))
                x_prior = x_prior + ...
                        obj.B * (u - obj.u0);                    
            end
            obj.xi = obj.Omega * x_prior;
        end
                
        function obj = update(obj, z)
            % Compute posterior information matrix
            invR = inv(obj.R);
            obj.Omega = obj.Omega + obj.H' * invR * obj.H;
            
            % Compute predicted measurement
            P = inv(obj.Omega);
            x = P * obj.xi;
            z_hat = obj.h0 + obj.H * (x - obj.x0);
            % Compute innovation
            innov = z - z_hat;
            % Compute posterior information vector
            obj.xi = obj.xi + obj.H' * invR * innov;                               
        end    
    end
end