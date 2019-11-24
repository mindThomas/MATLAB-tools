classdef BootstrapSISFilter
    % Sequential importance Sampling or Sequential Monte Carlo particle
    % filter implementation
    % The basis of this is the Sequential Importance Sampling algorithm
    % See "7.3.1 Sequential Importance Sampling (IS)" from Course ChM015x
	properties %(SetAccess = private)
        particles % state vectors of all particles (each column)
        weights
    end          
    properties %(Access = private)
        propagation_proposal_distribution % conditional distribution object from which samples can be drawn with .draw(x_prev)
        propagation_pdf % conditional PDF of the propagation distribution, p(x[k] | x[k-1]), on the form @(x, x_prev) ...
        measurement_pdf % conditional PDF of the measurement distribution, p(z[k] | x[k]), on the form @(z, x) ...        
    end  
    methods
        function obj = BootstrapSISFilter(propagation_proposal_distribution, propagation_pdf, measurement_pdf)                        
            obj.propagation_proposal_distribution = propagation_proposal_distribution;
            obj.propagation_pdf = propagation_pdf;
            obj.measurement_pdf = measurement_pdf;
        end        
        
        function obj = init(obj, n_particles, initial_distribution)
            % Initialize particle states by drawing an initial distribution and set a uniform weight            
            obj.particles = zeros(length(x_min), n_particles);
            obj.weights = 1/n_particles * ones(1, n_particles);
            for (j = 1:size(obj.particles, 2))
                obj.particles(:,j) = initial_distribution.draw();
            end
        end 
        
        function obj = filter(obj, z)
            for (j = 1:size(obj.particles, 2))
                % Propagate particles by drawing a new state from the
                % proposal distribution
                x_prev = obj.particles(:,j);
                x = obj.propagation_proposal_distribution.draw(x_prev);
                
                % Update the weight of the particle according to the
                % measurement likelihood multiplied with the proposal to
                % target distribution ratio
                likelihood = obj.measurement_pdf(z, x);
                proposal_ratio = obj.propagation_pdf(x, x_prev) / obj.propagation_proposal_distribution.pdf(x, x_prev);
                
                obj.particles(:,j) = x;
                obj.weights(j) = obj.weights(j) * likelihood * proposal_ratio;
            end             
            
            % Normalize weights
            obj.weights = obj.weights / sum(obj.weights);
        end               

        % Only for testing
        function obj = propagate(obj)
            for (j = 1:size(obj.particles, 2))                
                obj.particles(:,j) = obj.propagation_proposal_distribution.draw(x_prev);                
            end
        end               
        
        function plot(obj)
            if (size(obj.particles, 1) == 1)
                stem(obj.particles, obj.weights);
            end
        end
        
        function x_hat = getMMSEestimate(obj)
            % MMSE = Weighted average
            x_hat = zeros(size(obj.particles, 1), 1);
            for (j = 1:size(obj.particles, 2))
                x_hat = x_hat + obj.weights(j) * obj.particles(:,j);                 
            end
        end
        
        function x_hat = getMAPestimate(obj)
            % MAP = Find the most propable particle (highest weight)
            [val,idx] = max(obj.weights);
            x_hat = obj.particles(:,idx);
        end
        
        function mean = getNonlinearExpectation(obj, g)
            % Compute the expectation of a non-linear function with the
            % input set to the distribution modelled by the particle filter
            % E[g(x)]
            % where x ~ PF
            % and the function is defined as @(x) ...            
            test_input = zeros(size(obj.particles, 1), 1);
            mean = zeros(size(g(test_input), 1), 1);            
            for (j = 1:size(obj.particles, 2))
                mean = mean + obj.weights(j) * g(obj.particles(:,j)); 
            end
        end
            
    end
end