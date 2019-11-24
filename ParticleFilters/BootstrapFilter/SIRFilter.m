classdef SIRFilter < BootstrapSISFilter
    % Sequential importance Resampling filter implementation
    % The problem with the SIS filter is the degeneracy of the particles
    % (all but only one or a few particles have negligible=almost zero weight)
    % since the weights are just going down and down each each time
    % instance, since each particle has to capture the probability of the
    % whole particle trajectory to be true.
    % The SIR filter fixes this by implementing a Resambling technique
    % -> After resampling all particles will have equal weight
    % See "7.4.1 Sequential Importance Resampling (SIR)" from Course ChM015x
    %
    % Note that Particle filters suffers from the curse of dimensionality
    % and are intractable at higher dimensions
	%properties %(SetAccess = private)
    %    particles % state vectors of all particles (each column)
    %    weights
    %end          
    %properties (Access = private)
    %    propagation_proposal_distribution % conditional distribution object from which samples can be drawn with .draw(x_prev)
    %    propagation_pdf % conditional PDF of the propagation distribution, p(x[k] | x[k-1]), on the form @(x, x_prev) ...
    %    measurement_pdf % conditional PDF of the measurement distribution, p(z[k] | x[k]), on the form @(z, x) ...        
    %end  
    methods
        function obj = SIRFilter(propagation_proposal_distribution, propagation_pdf, measurement_pdf)                        
            obj = obj@BootstrapSISFilter(propagation_proposal_distribution, propagation_pdf, measurement_pdf);
        end        
                
        function obj = filter(obj, z)
            % Propagate and update using SIS Filter
            obj = filter@BootstrapSISFilter(obj, z);
            
            % The SIR filter is different in the way that the particles are
            % resampled after updating/recomputing the weights
            N = size(obj.particles, 2);
            if (obj.getEffectiveNumberOfParticles() < N/4)
                obj = obj.resample();
            end
        end
        
        function Neff = getEffectiveNumberOfParticles(obj)
            Neff = 1 / sum(obj.weights.^2);
            % Consider to resample when this number goes below N/4
        end
        
    end
    
    methods (Access = private)
        function obj = resample(obj)
            % Draw N samples with replacement from the current set of
            % particles with the probability of drawing a particular
            % particle defined according to its' weight.
            % We draw the samples by drawing a number from a uniform
            % distribution between 0 to 1 and then using a cumulated sum
            % lookup table to find the index of the according particle.
            new_particles = zeros(size(obj.particles));
                        
            resampling_wheel = cumsum(obj.weights);
            uniform_samples = rand(1, size(obj.particles, 2));
            for (i = 1:size(obj.particles, 2))                
                % Perform resampling wheel and look up in table to find
                % matching index
                sample = uniform_samples(i);
                old_particle_idx = length(resampling_wheel);
                for (j = 2:length(resampling_wheel))
                    if (resampling_wheel(j) > sample)
                        old_particle_idx = j;
                        break;
                    end
                end
                new_particles(:,i) = obj.particles(:,old_particle_idx);
            end
            obj.particles = new_particles;
            
            % Set all the weights and uniformly
            obj.weights = ones(size(obj.weights)) / size(obj.particles, 2);
        end
    end
end