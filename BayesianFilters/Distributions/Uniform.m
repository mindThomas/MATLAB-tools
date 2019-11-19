classdef Uniform
	properties (SetAccess = private)
        x_min
        x_max
    end    
    properties (Access = private)
        n
        scaleFactor  % used when computing products for handling un-normalized distributions
    end
    methods
        function obj = Uniform(x_min, x_max)
            if (size(x_min, 1) ~= size(x_max, 1))
                error('Inconsistent vector size');
            end
            
            % Uncorrelated Uniform distribution            
            obj.x_min = x_min;
            obj.x_max = x_max;
            obj.n = size(x_min, 1);  
            obj.scaleFactor = 1; % normalized
        end
        
        function x = draw(obj)
            % Draw a realization from the distribution
            x = (obj.x_max - obj.x_min) .* rand(obj.n,1) + obj.x_min;
        end
        
        function probability = pdf(obj, x)          
            % Compute the Gaussian PDF at a given point
            if (size(x, 1) ~= obj.n)
                error('Incorrect size')
            end
            
            prob = 1;
            for (i = 1:obj.n)
                prob = prob / (obj.x_max(i) - obj.x_min(i));
            end            
            
            probability = prob * ones(1, size(x, 2));
            for (j = 1:size(x, 2))
                for (i = 1:obj.n)
                    if ((x(i,j) < obj.x_min(i)) || (x(i,j) > obj.x_max(i)))
                        probability(1,j) = 0;
                    end
                end            
            end
            
            % Apply scale factor if distribution is not normalized
            probability = obj.scaleFactor * probability;
        end                
        
        function plot(obj, varargin)
            if (nargin == 2)
                fig = varargin{1};
                figure(fig);
            end

            if (obj.n == 1)
                % univariate distribution                
                d = obj.x_max - obj.x_min;
                x = [obj.x_min-d/10, obj.x_min, obj.x_max, obj.x_max+d/10];
                p = obj.pdf(obj.x_min) * [0, 1, 0, 0];                
                stairs(x, p);
                
            elseif (obj.n == 2)
                % bivariate distribution                
                d = obj.x_max - obj.x_min;
                xmin = obj.x_min(1) - d(1)/10;
                xmax = obj.x_max(1) + d(1)/10;
                ymin = obj.x_min(2) - d(2)/10;
                ymax = obj.x_max(2) + d(2)/10;
                steps = 50;
                xres = (xmax - xmin) / steps;
                yres = (ymax - ymin) / steps;
                x = (xmin:xres:xmax);
                y = (ymin:yres:ymax);

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
            end
        end        
    end    
    
    methods 
        function obj = normalize(obj)
            % Normalize distribution
            obj.scaleFactor = 1;
        end
        
        function joint = join(obj, g, varargin)
           % Join two Uniform distributions independently
           x_min2 = [obj.x_min; g.x_min];
           x_max2 = [obj.x_max; g.x_max];
           joint = Uniform(x_min2, x_max2);
        end                        
        
        function obj = marginalize(obj, mIdx)
            % Marginalize out the variables listed in the mIdx scalar or vector
            % from the variable vector: x = [a,b,c,d,...,r]
            % p(a,...,q) = int_r p(a,...,r) dr
            % For Multivariate Gaussians there is however an easy principle
            % that we can just drop the variables that we want to marginalize out
            % https://www.wikiwand.com/en/Multivariate_normal_distribution#/Marginal_distributions
            x_min2 = obj.x_min;
            x_min2(mIdx) = [];
            x_max2 = obj.x_max;
            x_max2(mIdx) = [];            
            obj = Uniform(x_min2, x_max2);
        end              
    end
end