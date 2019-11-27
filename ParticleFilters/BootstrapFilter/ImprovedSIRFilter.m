classdef ImprovedSIRFilter < SIRFilter
    % Sequential importance Resampling filter implementation
    % The problem with the SIS filter is the degeneracy of the particles
    % (all but only one or a few particles have negligible=almost zero weight)
    % since the weights are just going down and down each each time
    % instance, since each particle has to capture the probability of the
    % whole particle trajectory to be true.
    % The SIR filter fixes this by implementing a Resampling technique
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
    properties (Access = private)
        % To help with particle deprivation, particles are sampled from a
        % uniform distribution with probability, p, instead of being
        % sampled from the motion model
        % Particle deprivation = When performing estimation in a
        % high-dimensional space there may be no particles in the vicinity
        % to the correct state. Particle deprivation occurs as
        % the result of random resampling; an unlucky series of random 
        % numbers can wipe out all particles near the true state. 
        % In practice, problems of this nature only tend to arise when M is small relative
        % to the space of all states with high likelihood
        uniform_sampling_distribution
        uniform_sampling_probability
    end
    methods
        function obj = ImprovedSIRFilter(propagation_proposal_distribution, propagation_pdf, measurement_pdf, uniform_sampling_distribution, uniform_sampling_probability)
            obj = obj@SIRFilter(propagation_proposal_distribution, propagation_pdf, measurement_pdf);
            obj.uniform_sampling_distribution = uniform_sampling_distribution; % use this for 
            obj.uniform_sampling_probability = uniform_sampling_probability;
        end        
        
        function obj = filter(obj, z)
            uniform_samples = rand(size(obj.particles, 2), 1);
            for (j = 1:size(obj.particles, 2))
                % Propagate particles by drawing a new state from the
                % proposal distribution
                % However to reduce particle deprivation, draw a sample
                % from a general (preferably uniform) distribution with
                % probability (uniform_sampling_probability) to add some
                % random particles over the operating space
                x_prev = obj.particles(:,j);
                if (uniform_samples(j) <= obj.uniform_sampling_probability)
                    x = obj.uniform_sampling_distribution.draw();
                else
                    x = obj.propagation_proposal_distribution.draw(x_prev);
                end
                
                % Update the weight of the particle according to the
                % measurement likelihood multiplied with the proposal to
                % target distribution ratio
                likelihood = obj.measurement_pdf(z, x);
                proposal_pdf = obj.uniform_sampling_probability * obj.uniform_sampling_distribution.pdf(x) + ...
                               (1-obj.uniform_sampling_probability) * obj.propagation_proposal_distribution.pdf(x, x_prev);
                proposal_ratio = obj.propagation_pdf(x, x_prev) / proposal_pdf;                                
                
                obj.particles(:,j) = x;
                obj.weights(j) = obj.weights(j) * likelihood * proposal_ratio;
            end             
            
            % Normalize weights
            obj.weights = obj.weights / sum(obj.weights);        
            
            % The SIR filter is different in the way that the particles are
            % resampled after updating/recomputing the weights
            %
            % There exist two major strategies for variance reduction.
            % First, one may reduce the frequency at which resampling takes place.
            %
            % Resampling too often increases the risk of losing diversity.
            % If one samples too infrequently, many samples might be wasted
            % in regions of low probability. A standard approach to
            % determining whether or not resampling should be performed
            % is to measure the variance of the importance weights. 
            % If all weights are identical, then the variance is zero and
            % no resampling should be performed. If, on the other hand, the
            % weights are concentrated on a small number of samples, then
            % the weight variance is high and resampling should be performed
            N = size(obj.particles, 2);
            if (obj.getEffectiveNumberOfParticles() < N/4)
                disp('Resampling');
                %obj = obj.resample();
                %obj = obj.low_variance_resample();                
                obj = obj.efficient_low_variance_resample();
            end
        end                
        
        function Neff = getEffectiveNumberOfParticles(obj)
            Neff = 1 / sum(obj.weights.^2);
            % Consider to resample when this number goes below N/4
        end
        
    end
    
    methods (Access = private)
        function obj = low_variance_resample(obj)
            % Draw N samples from the current set as a sequential stochastic
            % process with the probability of drawing a particular particle
            % defined according to its' weight. Note that this is different
            % from the normal resampling strategy where samples are drawn
            % independently.
            %
            % Note also that this type of sampling does not require the
            % weights to be normalized (sum to one)
            %
            % The resampling step is a probabilistic implementation of the
            % Darwinian idea of survival of the fittest: It refocuses the
            % particle set to regions in state space with high posterior
            % probability
            %
            % We draw the samples by drawing ONE random number, r, from a
            % uniform distribution between 0 to 1/N and then using a cumulated sum
            % lookup table to find the index of the according particle.
            new_particles = zeros(size(obj.particles));
                        
            resampling_wheel = cumsum(obj.weights);
            N = size(obj.particles, 2);
            r = rand(1, 1) / N;  
            u = r;
            % OBS. That the below can be optimized (made more efficiently)
            % according to the algorithm "Algorithm Low variance sampler"
            % in Probabilistic Robotics by Sebastian Thrun
            for (i = 1:size(obj.particles, 2))                
                % Perform resampling wheel and look up in table to find
                % matching index starting from the selected random weight,
                % r, and moving 1/N at each time                                
                old_particle_idx = length(resampling_wheel);
                for (j = 2:length(resampling_wheel))
                    if (resampling_wheel(j) > u)
                        old_particle_idx = j;
                        break;
                    end
                end
                new_particles(:,i) = obj.particles(:,old_particle_idx);
                u = u + 1/N;
            end
            obj.particles = new_particles;
            
            % Set all the weights and uniformly
            obj.weights = ones(size(obj.weights)) / size(obj.particles, 2);
        end
        
        function obj = efficient_low_variance_resample(obj)
            % Draw N samples from the current set as a sequential stochastic
            % process with the probability of drawing a particular particle
            % defined according to its' weight. Note that this is different
            % from the normal resampling strategy where samples are drawn
            % independently.
            %
            % Note also that this type of sampling does not require the
            % weights to be normalized (sum to one)
            %
            % The resampling step is a probabilistic implementation of the
            % Darwinian idea of survival of the fittest: It refocuses the
            % particle set to regions in state space with high posterior
            % probability
            %
            % We draw the samples by drawing ONE random number, r, from a
            % uniform distribution between 0 to 1/N and then using a cumulated sum
            % lookup table to find the index of the according particle.
            new_particles = zeros(size(obj.particles));                                    
            N = size(obj.particles, 2);
            r = rand(1, 1) / N;  % draw random number between 0 to 1/N 
            u = r;
            accumulated_weight = obj.weights(:,1);
            particle_index = 1;            
            for (i = 1:size(obj.particles, 2))                
                % Perform resampling wheel and look up in table to find
                % matching index starting from the selected random weight,
                % r, and moving 1/N at each time                                
                while (u > accumulated_weight)
                    particle_index = particle_index + 1;
                    accumulated_weight = accumulated_weight + obj.weights(:,particle_index);
                end
                new_particles(:,i) = obj.particles(:,particle_index);
                u = u + 1/N;
            end
            obj.particles = new_particles;
            
            % Set all the weights and uniformly
            % Note that this step can be omitted when using this type of
            % sampling
            obj.weights = ones(size(obj.weights)) / size(obj.particles, 2);
        end          
    end
end