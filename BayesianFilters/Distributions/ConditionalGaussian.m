classdef ConditionalGaussian
	properties (SetAccess = private)
        Mu % mean function on the form: Mu = @(y) ...
        Cov % covariance matrix function: Cov = @(y) ...
    end    
    properties (Access = private)
        nx % number of elements in Gaussian vector
        ny % number of conditional variables     
        scaleFactor  % used when computing products for handling un-normalized distributions
    end
    methods
        function obj = ConditionalGaussian(mean, covariance, ny)  
            % Handle conditional Gaussian distributions of the form:
            % p(x|y)
            % where y is not only just one single value (e.g. y=1)
            % but the object instead captures the functional properties of
            % the conditional Gaussian
            
            % Compute test mean and covariance for checking size
            y = zeros(ny, 1);
            
            mean_test = mean(y);
            covariance_test = covariance(y);
            
            if (size(mean_test, 1) ~= size(covariance_test, 1) || size(mean_test, 1) ~= size(covariance_test, 2))
                error('Inconsistent mean and covariance size');
            end

            obj.Mu = mean;
            obj.Cov = covariance;
            obj.nx = size(mean_test, 1);
            obj.ny = ny;
            obj.scaleFactor = 1; % normalized
        end
        
        function x = draw(obj, y)
            % Draw a realization from the distribution given the
            % conditional, y
            mu = obj.Mu(y)
            C = obj.Cov(y);
            sqrtCov = chol(C, 'lower');           
            x = sqrtCov * randn(obj.n,1) + mu;
        end
        
        function probability = pdf(obj, x, y)          
            % Compute the Gaussian PDF at a given point given the
            % conditional, y
            mu = obj.Mu(y);
            C = obj.Cov(y);
            invCov = inv(C);
            detCov = det(C);
            if (size(x, 1) ~= obj.nx)
                error('Incorrect size')
            end
            if (size(x, 2) == 1)
                probability = 1/sqrt((2*pi)^obj.n * detCov) * exp(-1/2 * (x-mu)' * invCov * (x-mu));
            else
                probability = zeros(1, size(x, 2));
                for (i = 1:size(x, 2))
                    probability(i) = 1/sqrt((2*pi)^obj.n * detCov) * exp(-1/2 * (x(:,i)-mu)' * invCov * (x(:,i)-mu));
                end
            end
            
            % Apply scale factor if distribution is not normalized
            probability = obj.scaleFactor * probability;
        end        
        
        function likelihood = logLikelihood(obj, x, y)
            % Given an observation, x, and the conditional, y, what is the
            % log-likelihood of this observation
            mu = obj.Mu(y);
            C = obj.Cov(y);
            invCov = inv(C);
            detCov = det(C);
            likelihood = -1/2 * ( log(detCov) + (x-mu)'*invCov*(x-mu) + obj.n*log(2*pi) );
        end            
    end
end