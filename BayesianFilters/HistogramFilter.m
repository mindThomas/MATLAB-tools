classdef HistogramFilter < Histogram
    % Implementation of a Discrete Bayes filter working on discretized/quantized Histograms
    % See Table 4.1 in "Probabilistic Robotics" by Sebastian Thrun
    properties (Access = private)
        f % conditional PDF of the propagation distribution, p(x[k] | x[k-1]), on the form @(x, x_prev) ...
        h % measurement function
    end
    methods
        function obj = HistogramFilter(x_min, x_max, num_bins, f, h, init_pdf)           
            f_test = f(x_min, x_max);
            if (size(f_test,1) ~= size(x_min,1))
                error('Histogram filtering does not yet support change of domain');
            end
            
            h_test = h(x_min);
            if (size(h_test,1) ~= 1)
                error('Invalid measurement model PDF');
            end            
            
            obj = obj@Histogram(x_min, x_max, num_bins, init_pdf);
            obj.f = f;
            obj.h = h;
        end
        
        function obj = predict(obj)
            % Propagate all the discrete states bgased on the conditional PDF
            %    f(x, x_prev) = p(x[k] | x[k-1])            
            p_prev = obj.Probabilities;
            
            % Go through all discrete states and compute new propagated
            % probabilities
            for (k = 1:size(obj.BinCenters, 2))
                x = obj.BinCenters(:,k);
                obj.Probabilities(k) = 0;
                for (i = 1:size(obj.BinCenters, 2))
                    x_prev = obj.BinCenters(:,i);
                    obj.Probabilities(k) = obj.Probabilities(k) + ...
                       obj.f(x, x_prev) * p_prev(i);
                end
            end
            
            obj = obj.normalize();
        end
        
        function obj = update(obj, z)
            % Compute the likelihood of the measurement for all the
            % discrete states and update the probability according to Bayes
            % rule
            for (k = 1:size(obj.BinCenters,2))
                x = obj.BinCenters(:,k);
                obj.Probabilities(k) = h(x) * obj.Probabilities(k);
            end
            
            % Normalize probabilities
            obj = obj.normalize();
        end
        
    end
end