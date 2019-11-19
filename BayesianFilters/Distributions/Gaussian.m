classdef Gaussian
	properties (SetAccess = private)
        Mu
        Cov % covariance matrix
    end    
    properties (Access = private)
        n
        deterministic
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
            if (max(any(covariance)) == 0)
                obj.deterministic = true; % deterministic only supported if all variables as part of vector is deterministic                
                obj.sqrtCov = zeros(size(obj.Cov));
            else
                obj.deterministic = false;            
                obj.invCov = inv(covariance);
                obj.detCov = det(covariance);
                obj.sqrtCov = chol(covariance, 'lower');   
            end
            obj.scaleFactor = 1; % normalized
        end
        
        function x = draw(obj)
            % Draw a realization from the distribution
            if (obj.deterministic)
                x = obj.Mu;                
            else
                x = obj.sqrtCov * randn(obj.n,1) + obj.Mu;
            end
        end
        
        function probability = pdf(obj, x)          
            % Compute the Gaussian PDF at a given point
            if (size(x, 1) ~= obj.n)
                error('Incorrect size')
            end
            if (size(x, 2) == 1)
                if (obj.deterministic)
                    probability = 0;
                    return;
                end
                probability = 1/sqrt((2*pi)^obj.n * obj.detCov) * exp(-1/2 * (x-obj.Mu)' * obj.invCov * (x-obj.Mu));
            else
                probability = zeros(1, size(x, 2));
                if (obj.deterministic)
                    return;
                end
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
            if (obj.deterministic)
                likelihood = 0;
                return;
            end
            likelihood = -1/2 * ( log(obj.detCov) + (x-obj.Mu)'*obj.invCov*(x-obj.Mu) + obj.n*log(2*pi) );
        end    
        
        function [R, l] = getCovarianceRotation(obj)
            % Perform eigenvector decomposition to get covariance semiaxes
            % and their length
            if (obj.deterministic)
                R = eye(obj.n);
                l = zeros(obj.n,1);
                return;
            end
            
            [V, D] = eig(obj.Cov);
            E = sqrt(diag(D)); % the eigenvalues == sigma^2
            
            R = V;
            l = E;
        end
        
        function plot(obj, varargin)
            if (nargin == 2)
                fig = varargin{1};
                figure(fig);
            end

            if (obj.n == 1)
                % univariate distribution                
                x_min = obj.Mu - 4 * obj.sqrtCov;
                x_max = obj.Mu + 4 * obj.sqrtCov;
                steps = 200;
                res = (x_max - x_min) / steps;
                x = (x_min:res:x_max);
                p = obj.pdf(x);
                plot(x, p);
                
            elseif (obj.n == 2)
                % bivariate distribution                
                [R, l] = obj.getCovarianceRotation();
                
                four_sigma_points = [
                    obj.Mu + 4 * R(:,1) * l(1), ...
                    obj.Mu - 4 * R(:,1) * l(1), ...
                    obj.Mu + 4 * R(:,2) * l(2), ...
                    obj.Mu - 4 * R(:,2) * l(2)
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
                [R, l] = obj.getCovarianceRotation();
                
                for (s = 1:3)                
                    % generate data for "unrotated" ellipsoid
                    [xc,yc,zc] = ellipsoid(0,0,0,s*l(1),s*l(2),s*l(3));
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
        
        function plotSigmaContour(obj, sigma, varargin)
            if (obj.n ~= 2)
                error('This function only works with bivariate distributions');
            end
            if (nargin == 3)
                color = varargin{1};
            else
                color = 'r-';
            end

            % bivariate distribution            
            [R, l] = obj.getCovarianceRotation();

            % Generate ellipse
            theta = 0 : 0.01 : (2*pi+0.01);
            p = sigma * l(1) * R(:,1) * cos(theta) + sigma * l(2) * R(:,2) * sin(theta) + obj.Mu;
            plot(p(1,:), p(2,:), color);
                        
        end    
    end    
    
    methods 
        function obj = normalize(obj)
            % Normalize distribution
            obj.scaleFactor = 1;
        end
        
        function joint = join(obj, g, varargin)
           % Join two Gaussian distributions
           % varargin{1} = Pxy
           % if not specified Pxy = 0 and the distributions will be joint
           % as being independent
           if (nargin == 1)
               Pxy = vargin{1};
           else
               Pxy = zeros(size(obj.Mu,1), size(g.Mu,1));
           end
           Pyx = Pxy';
           Mu2 = [obj.Mu; g.Mu];
           Cov2 = [obj.Cov, Pxy;
                   Pyx, g.Cov];
           joint = Gaussian(Mu2, Cov2);
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
        
        function obj = conditional(obj, cvIdx, varargin)
            if (obj.deterministic)
                error('Conditionals are not possible with deterministic variables');
            end
            if (nargin == 3)                
                obj = conditional_evaluate(obj, cvIdx, varargin{1});
            else
                obj = conditional_generic(obj, cvIdx);
            end
        end                
            
        function obj2 = product(obj, g)
            % Compute the product between current Gaussian and g            
            % http://www.tina-vision.net/docs/memos/2003-003.pdf
            if (size(obj.Cov) ~= size(g.Cov))
                error('Distributions has to be equal in size');
            end
               
            if (obj.deterministic)
                obj2 = g;
                obj2.Mu = obj2.Mu * obj.Mu;
            elseif (g.deterministic)
                obj2 = obj;
                obj2.Mu = obj2.Mu * g.Mu;
            else
                invCov2 = obj.invCov + g.invCov;
                Cov2 = inv(invCov2);
                Mu2 = Cov2 * ( obj.invCov * obj.Mu + g.invCov * g.Mu);
                obj2 = Gaussian(Mu2, Cov2);
                obj2.scaleFactor = obj.pdf(Mu2) * g.pdf(Mu2) / obj2.pdf(Mu2);
            end
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
        
        function obj = add(obj, g)
            % Add two Gaussian random variables
            if (obj.n ~= g.n)
                error('Distributions not of the same size');
            end
            Mu2 = obj.Mu + g.Mu;
            Cov2 = obj.Cov + g.Cov;
            obj = Gaussian(Mu2, Cov2);
        end
        
        function obj2 = transform(obj, A, b)
            % Apply a linear transform to the random variables captured by
            % the distribution
            % Let the current distribution model p(x) = p(a,b,c,...)
            % This function transforms the distribution according to:
            % Z = A*x + b
            if (size(A,2) ~= size(obj.Mu,1))
                error('Invalid size of A');
            end            
            if (isnumeric(b))
                % Apply transformation with deterministic offset/constant, b
                if (size(b,1) == 0)
                    b = zeros(size(A,1),1);
                elseif (size(b,1) ~= size(A,1))
                    error('Invalid size of b');
                end               
                Mu2 = A*obj.Mu + b;
                % Cov(Z) = E[Z*Z'] - E[Z]*E[Z]'
                % Cov(Z) = E[(A*x + b)*(A*x + b)'] - (A*E[x]+b)*(A*E[x]+b)'
                % Cov(Z) = E[(A*x + b)*(x'*A' + b')] - (A*E[x]+b)*(E[x]'*A'+b')
                % Cov(Z) = E[(A*x + b)*(x'*A' + b')] - A*E[x]*E[x]'*A' - A*E[x]*b' - b*E[x]'*A' - b*b'
                % Cov(Z) = E[A*x*x'*A'] + E[b*x'*A'] + E[A*x*b'] + E[b*b'] - A*E[x]*E[x]'*A' - A*E[x]*b' - b*E[x]'*A' - b*b'
                % Cov(Z) = E[A*x*x'*A'] + b*E[x]'*A' + A*E[x]*b' + b*b' - A*E[x]*E[x]'*A' - A*E[x]*b' - b*E[x]'*A' - b*b'
                % Cov(Z) = A*E[x*x']*A' - A*E[x]*E[x]'*A'
                % Cov(Z) = A*(Cov(x)+E[x]*E[x]')*A' - A*E[x]*E[x]'*A'
                % Cov(Z) = A*Cov(x)*A'
                Cov2 = A * obj.Cov * A';
                obj2 = Gaussian(Mu2, Cov2);
            elseif (strcmp(class(b), 'Gaussian'))
                % Apply transformation with stochastic offset, b
                if (size(b.Mu,1) ~= size(A,1))
                    error('Invalid size of b');
                end                                
                Mu2 = A*obj.Mu + b.Mu;
                % Cov(Z) = E[Z*Z'] - E[Z]*E[Z]'
                % Cov(Z) = E[(A*x + b)*(A*x + b)'] - (A*E[x]+E[b])*(A*E[x]+E[b])'
                % Cov(Z) = E[(A*x + b)*(x'*A' + b')] - (A*E[x]+E[b])*(E[x]'*A'+E[b]')
                % Cov(Z) = E[(A*x + b)*(x'*A' + b')] - A*E[x]*E[x]'*A' - A*E[x]*E[b]' - E[b]*E[x]'*A' - E[b]*E[b]'
                % Cov(Z) = E[A*x*x'*A'] + E[b*x'*A'] + E[A*x*b'] + E[b*b'] - A*E[x]*E[x]'*A' - A*E[x]*E[b]' - E[b]*E[x]'*A' - E[b]*E[b]'
                % Cov(Z) = A*E[x*x']*A' + E[b*x']*A' + A*E[x*b'] + E[b*b'] - A*E[x]*E[x]'*A' - A*E[x]*E[b]' - E[b]*E[x]'*A' - E[b]*E[b]'
                % Assuming x and b to be independent
                % Cov(Z) = A*E[x*x']*A' + E[b]*E[x']*A' + A*E[x]*E[b'] + E[b*b'] - A*E[x]*E[x]'*A' - A*E[x]*E[b]' - E[b]*E[x]'*A' - E[b]*E[b]'
                % Cov(Z) = A*E[x*x']*A' + E[b*b'] - A*E[x]*E[x]'*A' - E[b]*E[b]'
                % Cov(Z) = A*(E[x*x'] - E[x]*E[x]')*A' + E[b*b'] - E[b]*E[b]'
                % Cov(Z) = A*Cov(x)*A' + Cov(b)
                Cov2 = A * obj.Cov * A' + b.Cov;
                obj2 = Gaussian(Mu2, Cov2);
            else
                error('Unknown type of b');
            end            
        end
        
        function joint = join_transform(obj, A, b)
            % Create a joint distribution, p(z), with a transformation of the
            % the random variable, x, defined by the current distribution.
            % That is, form a new random vector:
            % z = [x, y]
            % where the new random variable, y, is a transformation as:
            % y = A*x + b
            % We thus have:
            % z = [x]  =  [I]       +   [0]
            %     [y]     [A] * x       [I] * b
            % Which can also be condensed as:
            % z = F*x + G*b
            %   where
            %   F = [I; A]
            %   G = [0; I]
            % The mean is computed as:
            % E[z] = [ E[x]     =  [   E[x]        ]
            %          E[y] ]      [ A*E[x] + E[b] ]
            % The covariance is computed as:
            % Cov(z) = F*Cov(x)*F' + G*Cov(b)*G'            
            if (size(A,2) ~= size(obj.Mu,1))
                error('Invalid size of A');
            end            
            if (isnumeric(b))
                % Apply transformation with deterministic offset/constant, b
                if (size(b,1) == 0)
                    b = zeros(size(A,1),1);
                elseif (size(b,1) ~= size(A,1))
                    error('Invalid size of b');
                end   
                Mu2 = [ obj.Mu
                        A*obj.Mu + b ];   
                F = [eye(size(obj.Mu,1)); A];                
                Cov2 = F*obj.Cov*F';
            elseif (strcmp(class(b), 'Gaussian'))
                % Apply transformation with stochastic offset, b
                if (size(b.Mu,1) ~= size(A,1))
                    error('Invalid size of b');
                end                                   
                Mu2 = [ obj.Mu
                        A*obj.Mu + b.Mu ];
                F = [eye(size(obj.Mu,1)); A];
                G = [zeros(size(obj.Mu,1), size(b.Mu,1)); eye(size(b.Mu,1))];
                Cov2 = F*obj.Cov*F' + G*b.Cov*G';                    
            else
                error('Unknown type of b');
            end     
            joint = Gaussian(Mu2, Cov2);
        end
    end
    
    methods (Access = private)        
        function obj = conditional_evaluate(obj, cvIdx, cvValue)
            % Extract conditional Gaussian distribution based on index
            % within the variable vector: x = [a,b,c,d,...,r]
            % p(a,...,q | r=value) = p(a,...,r) / p(r)
            % https://stats.stackexchange.com/questions/30588/deriving-the-conditional-distributions-of-a-multivariate-normal-distribution
            % See also http://athenasc.com/Bivariate-Normal.pdf which
            % illustrates the bivariate case
            % See also "4.2.1 A Bayesian derivation of the Kalman filter"
            % of course ChM015x
            
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
        
        function obj = conditional_generic(obj, cvIdx)
            % Construct the dynamic/generic ConditionalGaussian object          
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

            % Given the conditional the mean and covariance is defined as:
            %Mu_new = Mu1 + Cov12 * inv(Cov22) * (cvValue - Mu2);
            %Cov_new = Cov11 - Cov12 * inv(Cov22) * Cov21;
            Mu_conditional = Mu1 - Cov12 * inv(Cov22) * Mu2;
            Mu_conditional_slope = Cov12 * inv(Cov22);
            Cov_conditional = Cov11 - Cov12 * inv(Cov22) * Cov21;
            
            obj = ConditionalGaussian(Mu_conditional, Mu_conditional_slope, Cov_conditional);
        end
    end
end