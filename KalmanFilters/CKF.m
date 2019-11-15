classdef CKF
	properties %(SetAccess = private)
        x        
    end    
	properties (SetAccess = private)        
        P        
        S
    end        
    properties (SetAccess = private)%(Access = private)
        f        
        h
        Q
        R           
        Q_sqrt
        R_sqrt 
    end   
    methods
        function obj = CKF(varargin)                        

        end        
     
        function obj = init_continuous(obj, f, Qc, ts, h, R, x_init, P_init)
            % Prediction model is on the form:
            %   dx = f(x, u, q)
            %   where q ~ N(0, Q)
            % which is discretized using Forward Euler into:
            %   x[k+1] = x[k] + ts*f(x[k], u[k], q[k])
            %   where q[k] ~ N(0, Q)
            % Measurement model is on the form:
            %   z = h(x, r)
            obj.f = @(x, u, q) ts*f(x, u, q); % Forward Euler discretization 
            obj.h = h;            
            obj.x = x_init;
            obj.P = P_init; 
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
            obj.x = x_init;
            obj.P = P_init;
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
            
            % Discrete prediction model on the form:
            %   x[k+1] = f(x[k], u[k], q[k])
            %   where q[k] ~ N(0, Q)  

            % Compute prior estimate by applying the unscented transform:            
            [obj.x, obj.P] = obj.CubatureTransform(obj.f, obj.x, u, obj.P, obj.Q);
        end
        
        function obj = update(obj, z)
            % Measurement model is on the form:
            %   z = h(x, r)          [ defined with h = @(x,r) ... ]
            %   where r ~ N(0, R)

            % Compute measurement estimate and innovation covariance estimate by applying the unscented transform:            
            [z_hat, obj.S, Pxz] = obj.CubatureTransform(obj.h, obj.x, [], obj.P, obj.Q);

            % Compute posterior
            obj.x = obj.x + Pxz * inv(obj.S) * (z - z_hat);
            obj.P = obj.P - Pxz * inv(obj.S) * Pxz';
        end
    end
    
    methods (Static, Access = private)              
        function [Mu_y, Cov_yy, Cov_xy] = CubatureTransform(f, Mu_x, u, P, Q)
            % y = f(x,u,q)
            %   where x ~ N(x, P)
            %   and   q ~ N(0, Q)
            % Or for measurement model
            % y = h(x,r)
            %   where x ~ N(x, P)
            %   and   r ~ N(0, R)   (where R is input into Q variable)
            % To support this method we augment the function to embed the
            % noise into the argument variable by extending it and
            % furthermore embeds the input value, u, into it
            % x_prime = [x; q]
            % y = f_prime(x_prime)
            %    where f_prime = f(x, u, q)
            nx = size(Mu_x,1);
            if (nargin(f) == 3)
                f_prime = @(x) f(x(1:nx), u, x(nx+1:end));
            else
                f_prime = @(x) f(x(1:nx), x(nx+1:end));
            end
            Mu_prime = [Mu_x; zeros(size(Q,1),1)];
            Cov_prime = [P, zeros(size(P,1), size(Q,2));
                         zeros(size(Q,1), size(P,2)), Q];

            % Process 1st and 2nd moment through a non-linear function using moment
            % matching
            % See "6.6.1 Sigma-point methods" of course ChM015x
            n = size(Mu_prime, 1);
            sqrtCov = chol(Cov_prime, 'lower');
            
            % Initialize the Sigma points according to the unscented transform                          
            % Cubature rule means that we only have Sigma points spread
            % at each +/- sigma locations. None in the mean.
            nSigmaPoints = 2*n;
            sigmaPoints = zeros(n, nSigmaPoints);
                                              
            for (i = 1:n)
                sigmaPoints(:,i) = Mu_prime + sqrt(n) * sqrtCov(:,i);
                sigmaPoints(:,i+n) = Mu_prime - sqrt(n) * sqrtCov(:,i);
            end
            weights = 1 / (2*n)   *  ones(1, nSigmaPoints);
            
            % Evaluate the non-linear function at the sigma points
            nf = size(f_prime(Mu_prime), 1);
            fSigmaPoints = zeros(nf, nSigmaPoints);
            for (i = 1:nSigmaPoints)
                fSigmaPoints(:, i) = f_prime(sigmaPoints(:, i));
            end
            
            % Compute the propagated mean and covariance
            Mu_y = fSigmaPoints * weights';
            Cov_yy = (fSigmaPoints - Mu_y) * diag(weights) * (fSigmaPoints - Mu_y)';
            Cov_x_prime_y = (sigmaPoints - Mu_prime) * diag(weights) * (fSigmaPoints - Mu_y)';
            Cov_xy = Cov_x_prime_y(1:nx,:);
        end
    end
end