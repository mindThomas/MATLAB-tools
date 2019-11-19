classdef Histogram
	properties (SetAccess = private)        
        BinWidth
        BinCenters
        Probabilities
    end 
    properties (Access = private)
        n
    end
    methods
        function obj = Histogram(x_min, x_max, num_bins, pdf)
            % pdf is a continuous PDF function on the form @(x) ...
            if (size(x_min, 1) ~= size(x_max, 1))
                error('Inconsistent vector sizes');
            end
            obj.n = size(x_min, 1);
            obj.BinWidth = (x_max - x_min) ./ num_bins;            
                        
            individual_bin_centers = cell(obj.n, 1);
            for (i = 1:obj.n)
                bin_edges = [x_min(i), ...
                                           x_min(i) + obj.BinWidth(i) * (1:(num_bins(i)-1)), ...
                                           x_max(i)];
                individual_bin_centers{i} = 1/2 * (bin_edges(1:end-1) + bin_edges(2:end));
            end
            
            obj.BinCenters = individual_bin_centers{1};
            for (i = 2:obj.n)
                obj.BinCenters = [repelem(obj.BinCenters, 1, size(individual_bin_centers{i},2));
                                  repmat(individual_bin_centers{i}, 1, size(obj.BinCenters,2))];
            end

            % Evaluate PDF at all bin centers
            obj.Probabilities = zeros(1, size(obj.BinCenters, 2));
            for (i = 1:size(obj.BinCenters, 2))
                obj.Probabilities(i) = pdf(obj.BinCenters(:,i));
            end

            % Normalize histogram
            obj.Probabilities = obj.Probabilities / sum(obj.Probabilities);
        end

        function plot(obj, varargin)
            if (nargin == 2)
                fig = varargin{1};
                figure(fig);
            end

            if (obj.n == 1)                   
                bar(obj.BinCenters, obj.Probabilities, 1);
            elseif (obj.n == 2)                                
                [X,Y] = meshgrid(unique(obj.BinCenters(1,:)), unique(obj.BinCenters(2,:)));                
                surf(X, Y, reshape(obj.Probabilities, size(X)));
                
                xlabel('X');
                ylabel('Y');
                zlabel('Probability');
            end
        end

        function obj = propagate(obj, f)
            f_test = f(obj.BinCenters(:,1));
            if (size(f_test,1) ~= size(obj.BinCenters,1))
                error('Propagation does not yet support change of domain');
            end

            % Non-linear propagation of the histogram bins, assuming the state 
            % to be deterministic when conditioned on each bin
            % f is a propagation function on the form @(x) ...
            propagated_centers = zeros(size(obj.BinCenters));
            original_probabilities = obj.Probabilities;
            for (i = 1:size(obj.BinCenters, 2))
                propagated_centers(:,i) = f(obj.BinCenters(:,i));
            end

            % Recompute bin probabilities
            obj.Probabilities = zeros(size(obj.Probabilities));
            for (i = 1:size(obj.BinCenters, 2))
                for (j = 1:size(propagated_centers, 2))
                    xmin = obj.BinCenters(:,i) - obj.BinWidth/2;
                    xmax = obj.BinCenters(:,i) + obj.BinWidth/2;
                    if (prod(propagated_centers(:,j) >= xmin) && prod(propagated_centers(:,j) <= xmax))
                        obj.Probabilities(i) = obj.Probabilities(i) + original_probabilities(j);
                    end
                end                
            end
        end
    end
end