classdef Gaussian
	properties (SetAccess = private)
        Mu
        Cov % covariance matrix
    end    
    properties (Access = private)
        n
        invCov
        detCov
        sqrtCov
        scaleFactor  % used when computing products for handling un-normalized distributions
    end
    methods
        function obj = Gaussian(mean, covariance)                        
            if (size(mean, 1) ~= size(covariance, 1) || size(mean, 1) ~= size(covariance, 2))
                error('Inconsistent mean and covariance size');
            end
            obj.Mu = mean;
            obj.Cov = covariance;
            obj.n = size(covariance, 1);        
            obj.invCov = inv(covariance);
            obj.detCov = det(covariance);
            obj.sqrtCov = chol(covariance, 'lower');   
            obj.scaleFactor = 1; % normalized
        end
        
        function x = draw(obj)
            % Draw a realization from the distribution
            x = obj.sqrtCov * randn(obj.n,1) + obj.Mu;
        end
        
        function probability = pdf(obj, x)          
            % Compute the Gaussian PDF at a given point
            if (size(x, 1) ~= obj.n)
                error('Incorrect size')
            end
            if (size(x, 2) == 1)
                probability = 1/sqrt((2*pi)^obj.n * obj.detCov) * exp(-1/2 * (x-obj.Mu)' * obj.invCov * (x-obj.Mu));
            else
                probability = zeros(1, size(x, 2));
                for (i = 1:size(x, 2))
                    probability(i) = 1/sqrt((2*pi)^obj.n * obj.detCov) * exp(-1/2 * (x(:,i)-obj.Mu)' * obj.invCov * (x(:,i)-obj.Mu));
                end
            end
            
            % Apply scale factor if distribution is not normalized
            probability = obj.scaleFactor * probability;
        end        
        
        function likelihood = logLikelihood(obj, x)
            % Given an observation, x, what is the log-likelihood of this
            % observation
            likelihood = -1/2 * ( log(obj.detCov) + (x-obj.Mu)'*obj.invCov*(x-obj.Mu) + obj.n*log(2*pi) );
        end    
        
        function [R, l] = getCovarianceRotation(obj)
            % Perform eigenvector decomposition to get covariance semiaxes
            % and their length
            [V, D] = eig(obj.Cov);
            E = diag(D);
            
            R = V;
            l = E;
        end
        
        function plot(obj, fig)
            if (obj.n == 1)
                % univariate distribution
                figure(fig);
                x_min = obj.Mu - 4 * obj.sqrtCov;
                x_max = obj.Mu + 4 * obj.sqrtCov;
                steps = 200;
                res = (x_max - x_min) / steps;
                x = (x_min:res:x_max);
                p = obj.pdf(x);
                plot(x, p);
                
            elseif (obj.n == 2)
                % bivariate distribution
                figure(fig);
                [R, l] = obj.getCovarianceRotation();
                
                four_sigma_points = [
                    obj.Mu + 4 * R(:,1) * sqrt(l(1)), ...
                    obj.Mu - 4 * R(:,1) * sqrt(l(1)), ...
                    obj.Mu + 4 * R(:,2) * sqrt(l(2)), ...
                    obj.Mu - 4 * R(:,2) * sqrt(l(2))
                ];
                
                x_min = min(four_sigma_points(1,:));
                x_max = max(four_sigma_points(1,:));
                y_min = min(four_sigma_points(2,:));
                y_max = max(four_sigma_points(2,:));
                steps = 50;
                x_res = (x_max - x_min) / steps;
                y_res = (y_max - y_min) / steps;
                x = (x_min:x_res:x_max);
                y = (y_min:y_res:y_max);

                p = zeros(length(x), length(y));
                for (i = 1:length(x))
                    v = [repmat(x(i), [1,length(y)])
                         y];
                    p(i, :) = obj.pdf(v);
                end
                
                [X, Y] = meshgrid(x, y);
                surf(X, Y, p');
                
                xlabel('X');
                ylabel('Y');
                zlabel('Probability');
                
                
            elseif (obj.n == 3)
                % multivariate 3D distribution
                figure(fig);
                [R, l] = obj.getCovarianceRotation();
                
                for (s = 1:3)                
                % generate data for "unrotated" ellipsoid
                [xc,yc,zc] = ellipsoid(0,0,0,s*sqrt(l(1)),s*sqrt(l(2)),s*sqrt(l(3)));
                % rotate data with rotation matrix R and center Mu
                a = kron(R(:,1),xc); b = kron(R(:,2),yc); c = kron(R(:,3),zc);
                data = a+b+c; n = size(data,2);
                x = data(1:n,:)+obj.Mu(1); y = data(n+1:2*n,:)+obj.Mu(2); z = data(2*n+1:end,:)+obj.Mu(3);
                if (s == 1)
                    mesh(x, y, z, 'edgecolor', 'b');   
                elseif (s == 2)
                    mesh(x, y, z, 'edgecolor', 'g');   
                elseif (s == 3)
                    mesh(x, y, z, 'edgecolor', 'r');   
                end
                hold on;
                end
                hold off;
                alpha(0.5);
                axis equal;                
            end
        end
    end    
    
    methods 
        function obj = normalize(obj)
            % Normalize distribution
            obj.scaleFactor = 1;
        end
        
        function obj = join_independent(obj, g)
           % Join two Gaussian distributions independently
           Mu2 = [obj.Mu; g.Mu];
           Cov2 = [obj.Sigma; zeros(size(obj.Cov,1), size(g.Cov,2));
                     zeros(size(g.Cov,1), size(obj.Cov,2)), g.Cov];
           obj = Gaussian(Mu2, Cov2);
        end                        
        
        function obj = marginalize(obj, mIdx)
            % Marginalize out the variables listed in the mIdx scalar or vector
            % from the variable vector: x = [a,b,c,d,...,r]
            % p(a,...,q) = int_r p(a,...,r) dr
            % For Multivariate Gaussians there is however an easy principle
            % that we can just drop the variables that we want to marginalize out
            % https://www.wikiwand.com/en/Multivariate_normal_distribution#/Marginal_distributions
            Mu2 = obj.Mu;
            Mu2(mIdx) = [];
            Cov2 = obj.Cov;
            Cov2(mIdx,:) = [];
            Cov2(:,mIdx) = [];
            obj = Gaussian(Mu2, Cov2);
        end        
                
        function obj = conditional(obj, cvIdx, cvValue)
            % Extract conditional Gaussian distribution based on index
            % within the variable vector: x = [a,b,c,d,...,r]
            % p(a,...,q | r=value) = p(a,...,r) / p(r)
            % https://stats.stackexchange.com/questions/30588/deriving-the-conditional-distributions-of-a-multivariate-normal-distribution
            % See also http://athenasc.com/Bivariate-Normal.pdf which
            % illustrates the bivariate case
            
            % For a conditional distribution we assume that we can perform
            % the given transformation to the random variable, given that
            % we now know the conditional variable. Lets simplify to:
            % p(x|y) = p(x,y) / p(y)
            % In that case we have our joint distribution p(x,y) but now we
            % know y. Since a joint Gaussian distribution describes the
            % relationship through the covariance (correlation) as in:
            % x = E[x] + sqrt(Cov(x,x)) * u + sqrt(Cov(x,y)) * v
            % y = E[y] + sqrt(Cov(y,x)) * u + sqrt(Cov(y,y)) * v
            % Knowing y will now allow us to turn the last equation around
            % to substitute v
            % v = inv(sqrt(Cov(y,y))) * (y - E[y] - sqrt(Cov(y,x)) * u)
            % Such that we get
            % Z = E[x] + sqrt(Cov(x,x)) * u + sqrt(Cov(x,y)) * inv(sqrt(Cov(y,y))) * (y - E[y] - sqrt(Cov(y,x)) * u)
            % And we isolate the deterministic parts
            % Z = E[x] + sqrt(Cov(x,y)) * inv(sqrt(Cov(y,y))) * (y - E[y]) + (sqrt(Cov(x,x)) - sqrt(Cov(x,y)) * inv(sqrt(Cov(y,y))) * sqrt(Cov(y,x))) * u
            % Now we isolate for the mean and variance
            % E[Z] = E[x] + sqrt(Cov(x,y)) * inv(sqrt(Cov(y,y))) * (y - E[y])
            % Cov(Z) = (sqrt(Cov(x,x)) - sqrt(Cov(x,y)) * inv(sqrt(Cov(y,y))) * sqrt(Cov(y,x))) * (sqrt(Cov(x,x)) - sqrt(Cov(x,y)) * inv(sqrt(Cov(y,y))) * sqrt(Cov(y,x)))'          
            %   - Cov(x,y) * inv(Cov(y,y)) * Cov(y,x)
            %   + Cov(x,x)
            
            Mu1 = obj.Mu;
            Mu1(cvIdx) = [];
            Mu2 = obj.Mu(cvIdx);
            Cov11 = obj.Cov;
            Cov11(cvIdx, :) = [];
            Cov11(:, cvIdx) = [];
            Cov12 = obj.Cov;
            Cov12(cvIdx,:) = [];
            Cov12 = Cov12(:,cvIdx);
            Cov22 = obj.Cov(cvIdx, cvIdx);
            Cov21 = obj.Cov;
            Cov21(:,cvIdx) = [];
            Cov21 = Cov21(cvIdx,:);            
            
            Mu_new = Mu1 + Cov12 * inv(Cov22) * (cvValue - Mu2);
            Cov_new = Cov11 - Cov12 * inv(Cov22) * Cov21;
            
            obj = Gaussian(Mu_new, Cov_new);
        end
        
        function obj2 = product(obj, g)
            % Compute the product between current Gaussian and g            
            % http://www.tina-vision.net/docs/memos/2003-003.pdf
            if (size(obj.Cov) ~= size(g.Cov))
                error('Distributions has to be equal in size');
            end
               
            invCov2 = obj.invCov + g.invCov;
            Cov2 = inv(invCov2);
            Mu2 = Cov2 * ( obj.invCov * obj.Mu + g.invCov * g.Mu);
            obj2 = Gaussian(Mu2, Cov2);
            obj2.scaleFactor = obj.pdf(Mu2) * g.pdf(Mu2) / obj2.pdf(Mu2);
        end
        
        function obj = convolution(obj, g)
            % Perform a convolution between current Gaussian and g
            % f conv g = int_0^x f(x - tau) * g(tau) dtau
            % The distribution of the sum of two random variables is
            % computed from the convolution of the two distributions
            % http://www.tina-vision.net/docs/memos/2003-003.pdf
            if (obj.n ~= 1 || g.n ~= 1)
                error('Currently the Convolution function only works for univariate Gaussians');
            end
            Mu2 = obj.Mu + g.Mu;
            Cov2 = obj.Cov + g.Cov;
            obj = Gaussian(Mu2, Cov2);
        end
    end
end