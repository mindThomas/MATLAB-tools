classdef BootstrapFilter
    % See "7.1.1 An introduction to particle filtering" from Course ChM015x
    % The simplest form of a Particle filter
	properties %(SetAccess = private)
        particles % state vectors of all particles (each column)
        weights
    end          
    properties (SetAccess = private)
        propagation_distribution    % conditional distribution    
        measurement_distribution    % conditional distribution
    end  
    methods
        function obj = BootstrapFilter(propagation_distribution, measurement_distribution)                        
            obj.propagation_distribution = propagation_distribution;
            obj.measurement_distribution = measurement_distribution;
        end        
        
        function obj = init(obj, n_particles, initial_distribution)
            % Initialize particle states by drawing an initial distribution and set a uniform weight            
            obj.particles = zeros(length(x_min), n_particles);
            obj.weights = 1/n_particles * ones(1, n_particles);
            for (j = 1:size(obj.particles, 2))
                obj.particles(:,j) = initial_distribution.draw();
            end
        end

        function obj = propagate(obj)
            for (j = 1:size(obj.particles, 2))
                obj.particles(:,j) = obj.propagation_distribution.draw(obj.particles(:,j));
            end
        end
        
        function obj = update(obj, z)
            for (j = 1:size(obj.particles, 2))
                likelihood = obj.measurement_distribution.pdf(z, obj.particles(:,j));
                obj.weights(j) = obj.weights(j) * likelihood                
            end
        end        
        
    end
end