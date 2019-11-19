classdef MomentPropagation 
    methods (Static)
        function [Mu2, Cov2] = Propagate_Linear(Mu, Cov, A, b, varargin)
            % y = A*x + b
        end
        
        function [Mu2, Cov2] = Propagate_Linearized(Mu, Cov, f)
            % y = f(x)
            % Which is assumed to be linearized to
            % y = f(Mu) + df/dx * (x - Mu)
            Mu2 = f(Mu);
            dfdx = numjacobian_real(f, Mu);
            Cov2 = dfdx * Cov * dfdx';
        end
        
        function [Mu2, Cov2, sigmaPoints, fSigmaPoints] = Propagate_Unscented(Mu, Cov, f, varargin)
            % y = f(x)
            % Process 1st and 2nd moment through a non-linear function using moment
            % matching
            % See "6.6.1 Sigma-point methods" of course ChM015x
            n = size(Mu, 1);
            sqrtCov = chol(Cov, 'lower');
            
            % Initialize the Sigma points according to the unscented transform
            if (nargin == 4)
                W0 = varargin{1};
            else
                W0 = 1 - n/3; % if no W0 is specified, set it as if the input distribution was Gaussian
            end       
            
            % Unscented transform Sigma points is one in the mean and one
            % at each +/- sigma locations
            nSigmaPoints = 2*n + 1;
            sigmaPoints = zeros(n, nSigmaPoints);
            
            sigmaPoints(:,1) = Mu;                        
            for (i = 1:n)
                sigmaPoints(:,1+i) = Mu + sqrt(n / (1-W0)) * sqrtCov(:,i);
                sigmaPoints(:,1+i+n) = Mu - sqrt(n / (1-W0)) * sqrtCov(:,i);
            end
            weights = (1 - W0) / (2*n)   *  ones(1, nSigmaPoints);
            weights(1) = W0;
            
            % Evaluate the non-linear function at the sigma points
            nf = size(f(Mu), 1);
            fSigmaPoints = zeros(nf, nSigmaPoints);
            for (i = 1:nSigmaPoints)
                fSigmaPoints(:, i) = f(sigmaPoints(:, i));
            end
            
            % Compute the propagated mean and covariance
            Mu2 = fSigmaPoints * weights';
            Cov2 = (fSigmaPoints - Mu2) * diag(weights) * (fSigmaPoints - Mu2)';
        end
        
        function [Mu2, Cov2, sigmaPoints, fSigmaPoints] = Propagate_Cubature(Mu, Cov, f)
            % y = f(x)
            % Process 1st and 2nd moment through a non-linear function using moment
            % matching
            % See "6.6.1 Sigma-point methods" of course ChM015x
            n = size(Mu, 1);
            sqrtCov = chol(Cov, 'lower');
            
            % Initialize the Sigma points according to the unscented transform                          
            % Cubature rule means that we only have Sigma points spread
            % at each +/- sigma locations. None in the mean.
            nSigmaPoints = 2*n;
            sigmaPoints = zeros(n, nSigmaPoints);
                        
            for (i = 1:n)
                sigmaPoints(:,i) = Mu + sqrt(n) * sqrtCov(:,i);
                sigmaPoints(:,i+n) = Mu - sqrt(n) * sqrtCov(:,i);
            end
            weights = 1 / (2*n)   *  ones(1, nSigmaPoints);
            
            % Evaluate the non-linear function at the sigma points
            nf = size(f(Mu), 1);
            fSigmaPoints = zeros(nf, nSigmaPoints);
            for (i = 1:nSigmaPoints)
                fSigmaPoints(:, i) = f(sigmaPoints(:, i));
            end
            
            % Compute the propagated mean and covariance
            Mu2 = fSigmaPoints * weights';
            Cov2 = (fSigmaPoints - Mu2) * diag(weights) * (fSigmaPoints - Mu2)';
        end        
    end
end   