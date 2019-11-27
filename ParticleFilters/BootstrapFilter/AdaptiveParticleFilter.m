classdef AdaptiveParticleFilter < SIRFilter
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
        % Adaptive Particle filter: KLD related variables
        % Computes the number of particles needed to guarantee that with probability 1-rho, the K-L distance between
        % the MLE and the true distribution is less than epsilon
        % See https://papers.nips.cc/paper/1998-kld-sampling-adaptive-particle-filters.pdf
        % See also code snippet: https://github.com/OpenSLAM-org/openslam_kld-sampling/blob/master/kld-sampling.cc
        max_error % upper bound on the K-L distance (distance between the MLE and the true distribution)
        z_value % quantile threshold
        min_particles
        % Based on a Histogram representation which quantizes the state
        % space into bins into which state sampling are counted
        BinWidth
        BinCenters
        BinOccupied        
    end
    methods

        function obj = initKLD(obj, x_min, x_max, num_bins, rho, epsilon, min_particles)
            % Computes the number of particles needed to guarantee that with probability 1-rho, the K-L distance between
            % the MLE and the true distribution is less than epsilon
            obj.max_error = epsilon; % upper bound on the K-L distance (distance between the MLE and the true distribution)
            obj.z_value = stdnormal_inv(1-rho); % compute the upper 1-rho quantile of the standard normal distribution
            obj.min_particles = min_particles;
            
            % Construct KLD bin centers and width based on num_bins vector
            % defining the number of bins along each state dimension
            obj.BinWidth = (x_max - x_min) ./ num_bins;            
                        
            n = size(x_min, 1);
            individual_bin_centers = cell(n, 1);
            for (i = 1:n)
                bin_edges = [x_min(i), ...
                                           x_min(i) + obj.BinWidth(i) * (1:(num_bins(i)-1)), ...
                                           x_max(i)];
                individual_bin_centers{i} = 1/2 * (bin_edges(1:end-1) + bin_edges(2:end));
            end
            
            obj.BinCenters = individual_bin_centers{1};
            for (i = 2:n)
                obj.BinCenters = [repelem(obj.BinCenters, 1, size(individual_bin_centers{i},2));
                                  repmat(individual_bin_centers{i}, 1, size(obj.BinCenters,2))];
            end
            
            % Set all bins to unopccupied
            obj.BinOccupied = true(1, size(obj.BinCenters, 2));
        end
                
        function obj = filter(obj, z)
            uniform_samples = rand(size(obj.particles, 2), 1);
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
            
            % In contrast to the variance minimizing approaches and also
            % approaches that minimize compute by determining when to
            % resample, the adaptive particle filter uses the KL bound to
            % determine how many number of particles to resample. Note however
            % that resampling always takes place in every single iteration
            obj = obj.KLD_resample();
        end                
        
        function Neff = getEffectiveNumberOfParticles(obj)
            Neff = 1 / sum(obj.weights.^2);
            % Consider to resample when this number goes below N/4
        end
        
    end
    
    methods (Access = private)
        function obj = KLD_resample(obj)
            % Draw an adaptive number of samples with replacement from the 
            % current set of particles with the probability of drawing a particular
            % particle defined according to its' weight.
            % The adaptive number of samples is computed according to the
            % KL bound
            %
            % The resampling step is a probabilistic implementation of the
            % Darwinian idea of survival of the fittest: It refocuses the
            % particle set to regions in state space with high posterior
            % probability
            %
            % We draw the samples by drawing a number from a uniform
            % distribution between 0 to 1 and then using a cumulated sum
            % lookup table to find the index of the according particle.
            obj = obj.KLD_Reset();
            new_particles = [];                        
            resampling_wheel = cumsum(obj.weights);                     
            n_particles = 0;
            while (n_particles < obj.KLD_GetKLbound())                
                % Perform resampling wheel and look up in table to find
                % matching index
                sample = rand(1, 1);
                old_particle_idx = length(resampling_wheel);
                for (j = 2:length(resampling_wheel))
                    if (resampling_wheel(j) > sample)
                        old_particle_idx = j;
                        break;
                    end
                end
                
                particle = obj.particles(:,old_particle_idx);
                obj = obj.KLD_Update(particle);
                new_particles(:,end+1) = particle;
                n_particles = n_particles + 1;
            end
            s = sprintf('Number of particles: %d\n', n_particles);
            disp(s);
            obj.particles = new_particles;
            
            % Set all the weights and uniformly
            obj.weights = ones(1, size(obj.particles, 2)) / size(obj.particles, 2);
        end
        
        function obj = KLD_Reset(obj)
            % Set all bins to unopccupied
            obj.BinOccupied = false(1, size(obj.BinCenters, 2));
        end
        
        function obj = KLD_Update(obj, sampled_particle)
            % Update the KLD grid/histogram by checking if the bin is empty
            % and if so, setting the bin to occupied
            % This algorithm can definitely be optimized
            for (i = 1:size(obj.BinCenters, 2))                
                xmin = obj.BinCenters(:,i) - obj.BinWidth/2;
                xmax = obj.BinCenters(:,i) + obj.BinWidth/2;
                % Determine if sampled particle belongs to bin
                if (sampled_particle >= xmin && sampled_particle <= xmax)
                    obj.BinOccupied(i) = true;
                    return;
                end
            end       
        end
        
        function n_samples = KLD_GetKLbound(obj)
            % Returns the number of optimal samples to generate to fulfil
            % the KL criteria
            n_occupied_bins = sum(obj.BinOccupied);
            if (n_occupied_bins < 2)
                n_samples = obj.min_particles;
                return;
            end
            % From https://papers.nips.cc/paper/1998-kld-sampling-adaptive-particle-filters.pdf
            n_samples = ceil( (n_occupied_bins-1)/(2*obj.max_error) * (1 - 2/(9.0*(n_occupied_bins-1)) + sqrt(2/(9.0*(n_occupied_bins-1))) * obj.z_value)^3 );
        end
               
    end
end