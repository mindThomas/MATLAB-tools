classdef AMCL
    % AMCL implements the  Adaptive Particle filter
    % See "AdaptiveParticleFilter.m"
    properties
        % Map
        map
        
        % Motion model
        dt
        motionModel
        
        % Particles
        particles % state vectors of all particles (each column)
        weights
        
        mean_pose
        
        % KLD related variables
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
        function obj = AMCL(map, dt, n_initial_particles, init_pose)       
            obj.map = map;
            
            % Initialize particles spread uniformly across the map and set a uniform weight            
            obj.particles = zeros(3, n_initial_particles);
            obj.weights = 1/n_initial_particles * ones(1, n_initial_particles);
            for (j = 1:size(obj.particles, 2))
                %obj.particles(1,j) = (map.pmax(1)-map.pmin(1))*rand(1,1) + map.pmin(1); % x
                %obj.particles(2,j) = (map.pmax(2)-map.pmin(2))*rand(1,1) + map.pmin(2); % y
                %obj.particles(3,j) = (2*pi)*rand(1,1) - pi; % psi (yaw)
                obj.particles(1,j) = 1.0*rand(1,1) - 0.5 + init_pose(1); % x
                obj.particles(2,j) = 1.0*rand(1,1) - 0.5 + init_pose(2); % y
                obj.particles(3,j) = 0.5*rand(1,1) - 0.25 + init_pose(3); % psi (yaw)                
            end
            
            x_min = [-0.5; -0.5; -pi/2];
            x_max = [0.5; 0.5; pi/2]; 
            num_bins = [20; 20; 10];
            
            % Computes the number of particles needed to guarantee that with probability 1-rho, the K-L distance between
            % the MLE and the true distribution is less than epsilon
            rho = 0.05; % 5% guaranteed
            epsilon = 0.05;            
            min_particles = 50;
            obj = obj.initKLD(x_min, x_max, num_bins, rho, epsilon, min_particles);
            
            obj.dt = dt;            
            [f, Fx, Fu, Fq] = CoordinatedTurnModel_Discrete_Thrun(dt);
            obj.motionModel = f;
            
            obj.mean_pose = mean(obj.particles,2);
        end
                
        function obj = filter(obj, v, omega, range_measurements, sensors)
            % scan contains a cell array of BeamRangeSensor objects            
            log_likelihood = zeros(size(obj.particles, 2), 1);
            for (j = 1:size(obj.particles, 2))
                % Propagate particles by drawing a new state from the
                % proposal distribution
                x_prev = obj.particles(:,j);
                x = obj.drawFromMotionModel(x_prev, v, omega);                
                
                % Update the weight of the particle according to the
                % measurement likelihood multiplied with the proposal to
                % target distribution ratio
                log_likelihood(j) = 0;
                for (i = 1:length(range_measurements))                    
                    log_likelihood(j) = log_likelihood(j) + log(sensors{i}.inverse_measurement_probability(range_measurements(i), x, obj.map));
                end                
                %proposal_ratio = obj.propagation_pdf(x, x_prev) / obj.propagation_proposal_distribution.pdf(x, x_prev);                                
                %proposal_ratio = 1; % set to one assuming our proposal and target distribution are the same, even though they are not!
                
                obj.particles(:,j) = x;
            end             
            
            % Update weights
            max_log_likelihood = max(log_likelihood)
            for (j = 1:size(obj.particles, 2))
                obj.weights(j) = obj.weights(j) * exp(log_likelihood(j) + max_log_likelihood);
            end
            
            % Normalize weights
            obj.weights = obj.weights / sum(obj.weights);        
            
            % In contrast to the variance minimizing approaches and also
            % approaches that minimize compute by determining when to
            % resample, the adaptive particle filter uses the KL bound to
            % determine how many number of particles to resample. Note however
            % that resampling always takes place in every single iteration
            obj = obj.KLD_resample();
            
            obj.mean_pose = mean(obj.particles,2);
        end   
        
        function obj = propagate(obj, v, omega)
            % scan contains a cell array of BeamRangeSensor objects            
            for (j = 1:size(obj.particles, 2))
                % Propagate particles by drawing a new state from the
                % proposal distribution                
                obj.particles(:,j) = obj.drawFromMotionModel(obj.particles(:,j), v, omega);                                
            end             
        end        
        
        function obj = update(obj, range_measurements, sensors)
            % scan contains a cell array of BeamRangeSensor objects            
            for (j = 1:size(obj.particles, 2))
                % Update the weight of the particle according to the
                % measurement likelihood multiplied with the proposal to
                % target distribution ratio
                % Using log-likelihood for numerical reasons
                % See https://www.cs.utexas.edu/~kuipers/handouts/S07/L5%20Markov%20localization.pdf
                log_likelihood(j) = 0;
                for (i = 1:length(range_measurements))                    
                    log_likelihood(j) = log_likelihood(j) + log(sensors{i}.inverse_measurement_probability(range_measurements(i), obj.particles(:,j), obj.map));
                end               
                %proposal_ratio = obj.propagation_pdf(x, x_prev) / obj.propagation_proposal_distribution.pdf(x, x_prev);                                
                %proposal_ratio = 1; % set to one assuming our proposal and target distribution are the same, even though they are not!
            end     
            
            % Update weights
            max_log_likelihood = max(log_likelihood);            
            for (j = 1:size(obj.particles, 2))
                obj.weights(j) = obj.weights(j) * exp(log_likelihood(j) + max_log_likelihood);
            end            
            
            % Normalize weights
            obj.weights = obj.weights / sum(obj.weights);        
            
            % In contrast to the variance minimizing approaches and also
            % approaches that minimize compute by determining when to
            % resample, the adaptive particle filter uses the KL bound to
            % determine how many number of particles to resample. Note however
            % that resampling always takes place in every single iteration
            obj = obj.KLD_resample();
            
            obj.mean_pose = mean(obj.particles,2);
        end           
        
    end
    
    methods (Access = private)
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
        
        function newPose = drawFromMotionModel(obj, pose, v, omega)
            % Implement bicycle model for propagation
            % v = translational velocity
            % rho = steering angle
            % omega = angular velocity            

            % Correlation parameters
            alpha1 = 0.5; % velocity contribution into velocity variance 
            alpha2 = 0; % angular velocity contribution into velocity variance 
            alpha3 = 1.5; % velocity contribution into angular velocity variance 
            alpha4 = 10; % angular velocity contribution into angular velocity variance 
            alpha5 = 0; % velocity contribution into angle pertubation variance 
            alpha6 = 10; % angular velocity contribution into angle pertubation variance 
            % Define pertubation variables
            v_var = (alpha1*v^2 + alpha2*omega^2);
            omega_var = (alpha3*v^2 + alpha4*omega^2);
            gamma_var = (alpha5*v^2 + alpha6*omega^2);            
            gamma = 0;
            
            % Add noise            
            v = v + sqrt(v_var)*randn(1,1);
            omega = omega + sqrt(omega_var)*randn(1,1);
            gamma = sqrt(gamma_var)*randn(1,1);            
            
            % Practical fix against zero angular velocity which will cause
            % division by zero in the motion model
            if (abs(omega) < 100*eps)
                omega = 100*eps;
            end
            
            % Construct motion model inputs
            u = [v; omega];
            q = gamma;
            
            newPose = obj.motionModel(pose, u, q);
        end    
        
        function Neff = getEffectiveNumberOfParticles(obj)
            Neff = 1 / sum(obj.weights.^2);
            % Consider to resample when this number goes below N/4
        end
        
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
                if (prod((sampled_particle-obj.mean_pose) >= xmin) && prod((sampled_particle-obj.mean_pose) <= xmax))
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