classdef ConditionalNonlinearWithAdditiveGaussian
    properties (SetAccess = private)
        fmean % mean function conditioned on the conditional variable in the form of E[x] = @(y) ...
        Cov
    end
    properties (Access = private)        
        nx
        invCov
        detCov
        sqrtCov    
    end
    methods
        function obj = ConditionalNonlinearWithAdditiveGaussian(fmean, covariance)  
            % Handle non-linear conditional distributions with additive
            % Gaussian noise on the form
            % p(x|y)
            % where the output variable x is defined as:
            % x = f(y) + q
            % and q ~ N(0, cov)        
            obj.fmean = fmean;
            obj.nx = size(covariance, 1);
            obj.Cov = covariance;
            obj.invCov = inv(covariance);
            obj.detCov = det(covariance);            
            obj.sqrtCov = chol(covariance, 'lower');           
        end
        
        function obj2 = given(obj, y)
            % Get the distribution given a certain value of the conditional
            % variable. In this case the resulting distribution becomes
            % Gaussian since the variables entering the non-linear part is given.
            mu = obj.fmean(y);
            obj2 = Gaussian(mu, obj.Cov);
        end
        
        function x = draw(obj, y)
            % Draw a realization from the distribution given the
            % conditional, y
            mu = obj.fmean(y);     
            x = obj.sqrtCov * randn(obj.nx,1) + mu;
        end
        
        function probability = pdf(obj, x, y)          
            % Compute the Gaussian PDF at a given point given the
            % conditional, y
            mu = obj.fmean(y);
            if (size(x, 1) ~= obj.nx)
                error('Incorrect size')
            end
            if (size(x, 2) == 1)
                probability = 1/sqrt((2*pi)^obj.nx * obj.detCov) * exp(-1/2 * (x-mu)' * obj.invCov * (x-mu));
            else
                probability = zeros(1, size(x, 2));
                for (i = 1:size(x, 2))
                    probability(i) = 1/sqrt((2*pi)^obj.nx * obj.detCov) * exp(-1/2 * (x(:,i)-mu)' * obj.invCov * (x(:,i)-mu));
                end
            end
        end        
        
        function likelihood = logLikelihood(obj, x, y)
            % Given an observation, x, and the conditional, y, what is the
            % log-likelihood of this observation
            mu = obj.fmean(y);
            likelihood = -1/2 * ( log(obj.detCov) + (x-mu)'*obj.invCov*(x-mu) + obj.nx*log(2*pi) );
        end       
        
        function plot(obj, varargin)
            if (nargin == 2)
                y = varargin{1};
            else
                y = 0;
            end
            obj.given(y).plot();
        end
    end
end