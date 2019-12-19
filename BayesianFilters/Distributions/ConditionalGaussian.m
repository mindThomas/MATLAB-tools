classdef ConditionalGaussian
	properties (SetAccess = private)
        Mu
        Cov
    end    
    properties (Access = private)
        nx % number of elements in Gaussian vector
        ny % number of conditional variables     
        MuSlope % slope of mean with regards to the conditional variables
        invCov
        detCov
        sqrtCov
        scaleFactor  % used when computing products for handling un-normalized distributions
    end
    methods
        function obj = ConditionalGaussian(mean, mean_conditional_slope, covariance)  
            % Handle conditional Gaussian distributions of the form:
            % p(x|y)
            % where y is not only just one single value (e.g. y=1)
            % but the object instead captures the functional properties of
            % the conditional Gaussian           
            obj.Mu = mean;
            obj.MuSlope = mean_conditional_slope;
            obj.Cov = covariance;
            obj.invCov = inv(covariance);
            obj.detCov = det(covariance);
            obj.sqrtCov = chol(covariance, 'lower');   
            obj.nx = size(obj.Mu,1);
            obj.ny = size(obj.MuSlope,2);
            obj.scaleFactor = 1; % normalized
        end
        
        function obj2 = given(obj, y)
            % Get the distribution given a certain value of the conditional
            % variable
            mu = obj.Mu + obj.MuSlope * y;
            obj2 = Gaussian(mu, obj.Cov);
        end
        
        function x = draw(obj, y)
            % Draw a realization from the distribution given the
            % conditional, y
            mu = obj.Mu + obj.MuSlope * y;     
            x = obj.sqrtCov * randn(obj.nx,1) + mu;
        end
        
        function probability = pdf(obj, x, y)          
            % Compute the Gaussian PDF at a given point given the
            % conditional, y
            mu = obj.Mu + obj.MuSlope * y;
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
            
            % Apply scale factor if distribution is not normalized
            probability = obj.scaleFactor * probability;
        end        
        
        function likelihood = logLikelihood(obj, x, y)
            % Given an observation, x, and the conditional, y, what is the
            % log-likelihood of this observation
            mu = obj.Mu + obj.MuSlope * y;
            likelihood = -1/2 * ( log(obj.detCov) + (x-mu)'*obj.invCov*(x-mu) + obj.nx*log(2*pi) );
        end       
        
        function joint = join(obj, y)
            % Multiply this conditional distribution, p(x|y), with the
            % distribution of the conditional variable, p(y), to get the
            % joint distribution:
            % p(x,y) = p(x|y) * p(y)
            % Since this conditional is modelled as a linear transform
            % (which is the case when dealing with Gaussian conditionals)
            % Lets denote the conditional variable as z such that we denote
            % p(z) = p(x|y)
            % Now given the linear relationship, this variable z is defined as
            % z = x + slope*y
            % Hence the joint distribution variables can be modelled as
            % g = [ x + slope*y ]   =   [I]      +  [slope]
            %     [ y           ]       [0] * x     [  I  ] * y
            % Which can also be condensed as:
            % g = A*x + B*b
            %   where
            %   A = [I; 0]
            %   B = [slope; I]
            % The mean is computed as:
            % E[g] = [ E[x + slope*y]     =  [ E[x] + slope*E[y]  ]
            %          E[y]           ]      [ E[y]               ]
            % The covariance is computed as:
            % Cov(g) = A*Cov(x)*A' + B*Cov(y)*B'            
            if (~strcmp(class(y), 'Gaussian'))
                error('Conditional distribution should also be a Gaussian');
            end
            
            A = [eye(obj.nx)
                 zeros(obj.ny, obj.nx)];
            B = [obj.MuSlope
                 eye(size(obj.MuSlope,2))];
                                                  
            Mu2 = [ obj.Mu + obj.MuSlope * y.Mu
                    y.Mu ];
            Cov2 = A*obj.Cov*A' + B*y.Cov*B';
   
            joint = Gaussian(Mu2, Cov2);            
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