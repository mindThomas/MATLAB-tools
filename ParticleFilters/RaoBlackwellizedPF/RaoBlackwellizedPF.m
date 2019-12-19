classdef RaoBlackwellizedPF
    % Rao-Blackwellized Particle filter implementation
    % Also known as the Marginalized Particle filter
    % The Rao-Blackwellized PF is based upon the assumption that the states
    % of the distribution to estimate can be partitioned into a set of
    % linear states and non-linear states. Being that the linear states
    % only acts linearly into both propagation and measurement models.
    % xi = [x, y]
    % where x is the linear states
    % and y is the non-linear states
    % The propagation model should thus be given on the form
    % x[k] = f(y[k-1]) + Fx(y[k-1]) * x[k-1] + qx[k-1]
    % y[k] = g(y[k-1]) + Gx(y[k-1]) * x[k-1] + qy[k-1]
    % z[k] = h(y[k]) + Hx(y[k]) * x[k] + r[k]
    % where qx, qy and r are random Gaussian noise
    % Note that the version implemented here uses the motion model
    % propagation distribution as the proposal distribution
    
    % See "7.6.1 Rao-Blackwellized particle filter" from Course ChM015x
    % See also http://user.it.uu.se/~thosc112/pubpdf/schongk2011.pdf
    % And http://user.it.uu.se/~thosc112/pubpdf/schongn2005.pdf
    %
    % See also Chapter 7.5 of "Bayesian Filtering and Smoothing" from
    % https://users.aalto.fi/~ssarkka/pub/cup_book_online_20131111.pdf
    %
    % See this paper for improvements to the Proposal distribution for the
    % RBPF: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.73.1897&rep=rep1&type=pdf
	properties %(SetAccess = private)
        particles % cell-array of particle objects (RaoBlackwellizedParticle)        
    end          
    properties %(Access = private)
        n_linear
        n_nonlinear
        % Propagation models
        f % linear state propagation mean contribution conditioned on the non-linear variable, on the form @(y) ...
        Fx % linear state propagation transition matrix wrt. the linear states, conditioned on the non-linear variable, on the form @(y) ...
        g % non-linear state propagation transition function, on the form @(y) ...
        Gx % non-linear state propagation transition matrix wrt. the linear states, conditioned on the non-linear variable, on the form @(y) ...
        Cov_qxy % joint covariance of the additive Gaussian noise        
        % Measurement model
        h % measurement model mean contribution conditioned on the non-linear variable, on the form @(y) ...
        Hx % linear measurement model Jacobian conditioned on the non-linear variable, on the form @(y) ...
        Cov_r
        measurement_nonlinear_pdf % conditional PDF of the measurement distribution, p(z[k] | y[k]), on the form @(z, y) ...
                                  % note that if the measurement model depends on the linear states, the model will have to be
                                  % marginalized to remove this dependency on x and thereby make a likehood function independent of x
    end  
    methods
        function obj = RaoBlackwellizedPF(n_linear, f, Fx, n_nonlinear, g, Gx, Cov_qxy, h, Hx, Cov_r, measurement_nonlinear_pdf)
            obj.n_linear = n_linear;            
            obj.f = f;
            obj.Fx = Fx;
            obj.n_nonlinear = n_nonlinear;
            obj.g = g;
            obj.Gx = Gx;
            obj.Cov_qxy = Cov_qxy;
            obj.h = h;
            obj.Hx = Hx;
            obj.Cov_r = Cov_r;            
            obj.measurement_nonlinear_pdf = measurement_nonlinear_pdf;
        end        
        
        function obj = init(obj, n_particles, initial_nonlinear_distribution, initial_linear_mean, initial_linear_covariance)
            if (size(initial_linear_mean, 1) ~=  obj.n_linear || size(initial_linear_mean, 1) ~= size(initial_linear_covariance, 1) || size(initial_linear_mean, 1) ~= size(initial_linear_covariance, 2))
                error('Inconsistent mean and covariance size');
            end
            % Initialize particle states by drawing an initial distribution and set a uniform weight            
            obj.particles = {};
            %obj.weights = 1/n_particles * ones(1, n_particles);
            weight = 1/n_particles;
            for (j = 1:n_particles)                
                obj.particles{j} = RaoBlackwellizedParticle(weight, initial_nonlinear_distribution.draw(), initial_linear_mean, initial_linear_covariance);                
            end
        end 
        
        function obj = filter(obj, z)
            weights_sum = 0;
            
            for (j = 1:length(obj.particles))
                % 4. Non-linear state Particle filter update
                % Use the measurement model likelihood function to update the
                % weight of each particle. Note that if the measurement model
                % depends on the linear state(s), these will have to be
                % marginalized out when calculating this likelihood.
                y = obj.particles{j}.nonlinear_state;
                likelihood = obj.measurement_nonlinear_pdf(z, y);              
                % This could also have been defined/computed different, see
                % e.g.: http://user.it.uu.se/~thosc112/pubpdf/schongn2005.pdf
                % With code here: http://user.it.uu.se/~thosc112/research/rao-blackwellized-particle.html
                obj.particles{j}.weight = obj.particles{j}.weight * likelihood;
                weights_sum = weights_sum + obj.particles{j}.weight;               
            end
            
            % Normalize weights
            for (j = 1:length(obj.particles))
                obj.particles{j}.weight = obj.particles{j}.weight / weights_sum;
            end                        
                
            % Check for resampling need            
            N = length(obj.particles);
            if (obj.getEffectiveNumberOfParticles() < N/4)
                disp('Resampling');
                obj = obj.low_variance_resample();
            end                       
                 
            for (j = 1:length(obj.particles))
                % 5. Linear state Kalman filter measurement update
                % Perform a measurement update of the linear states using the
                % measurement model where the non-linear states are given by
                % the particle state (viewed as a known constant/deterministic).                
                h = obj.h(obj.particles{j}.nonlinear_state);
                Hx = obj.Hx(obj.particles{j}.nonlinear_state);     
                if (any(Hx))
                    obj.particles{j} = obj.particles{j}.linear_update(h, Hx, obj.Cov_r, z);                            
                end
                
                % 1. Non-linear state Particle filter prediction
                % Draw samples from the non-linear proposal density/distribution
                % Note that when drawing samples from the non-linear
                % distribution the uncertainty of the linear state has to be
                % considered by drawing a sample of the linear state and use
                % this when drawing a sample from the conditional distribution
                % for the non-linear propgation model            
                y_prev = obj.particles{j}.nonlinear_state;
                                
                % Compute resulting noise covariance of particle by
                % inflating the noise with the uncertainty in the linear state                
                cov_nonlinear = obj.Gx(y_prev) * obj.particles{j}.P * obj.Gx(y_prev)' + obj.Cov_qxy(obj.n_linear+1:end,obj.n_linear+1:end);                                
                % Draw noise sample
                qy = chol(cov_nonlinear, 'lower') * randn(obj.n_nonlinear,1);            
                %qy = chol(obj.Cov_qxy(obj.n_linear+1:end,obj.n_linear+1:end), 'lower') * randn(obj.n_nonlinear,1);                            
                % Compute new non-linear state
                y = obj.g(y_prev) + obj.Gx(y_prev) * obj.particles{j}.x + qy;
                obj.particles{j}.nonlinear_state = y;
                obj.particles{j}.prev_nonlinear_state = y_prev;                        

                                
                f = obj.f(obj.particles{j}.prev_nonlinear_state);
                Fx = obj.Fx(obj.particles{j}.prev_nonlinear_state);
                g = obj.g(obj.particles{j}.prev_nonlinear_state);
                Gx = obj.Gx(obj.particles{j}.prev_nonlinear_state);
                cov_qx = obj.Cov_qxy(1:obj.n_linear,1:obj.n_linear);
                cov_qy = obj.Cov_qxy(obj.n_linear+1:end,obj.n_linear+1:end);
                cov_qxy = obj.Cov_qxy(1:obj.n_linear,obj.n_linear+1:end);
                obj.particles{j} = obj.particles{j}.dynamics_update(obj.particles{j}.nonlinear_state, f, g, Fx, Gx, cov_qx, cov_qy, cov_qxy);
                
                
%                 % 3. Linear state Kalman filter prediction
%                 % Apply the linear motion model to predict/propagate the linear
%                 % state and covariance                            
%                 cov_qx = obj.Cov_qxy(1:obj.n_linear,1:obj.n_linear);
%                 cov_qy = obj.Cov_qxy(obj.n_linear+1:end,obj.n_linear+1:end);
%                 cov_qx_qy = obj.Cov_qxy(1:obj.n_linear,obj.n_linear+1:end);
%                 z_dyn = obj.particles{j}.nonlinear_state - obj.g(obj.particles{j}.prev_nonlinear_state);
%                 bias = obj.f(obj.particles{j}.nonlinear_state) + cov_qx_qy'*inv(cov_qy)*z_dyn; % handle the correlation between qx and qy
%                 A = obj.Fx(obj.particles{j}.nonlinear_state) - cov_qx_qy'*inv(cov_qy)*obj.Gx(obj.particles{j}.prev_nonlinear_state); % handle the correlation between qx and qy
%                 Q = cov_qx - cov_qx_qy'*inv(cov_qy)*cov_qx_qy;  % handle the correlation between qx and qy
%                 obj.particles{j} = obj.particles{j}.linear_predict(bias, A, Q);                     
%                 
%                 % 2. Linear state Kalman filter update with change in the non-linear state
%                 % The sample that was just drawn for the non-linear state is
%                 % used in an update step of the Kalman filter for the linear
%                 % state. This thus requires that an dynamics-based measurement
%                 % model is constructed for the linear state, relating a
%                 % previous non-linear state to the current new non-linear state
%                 % Note that this can be seen as using the motion model of the
%                 % non-linear state as a measurement model for the linear state.
%                 % y[k] = g(y[k-1]) + Gx(y[k-1]) * x[k-1] + qy[k-1]
%                 % y[k] - g(y[k-1]) = Gx(y[k-1]) * x[k-1] + qy[k-1]            
%                 z_dyn = obj.particles{j}.nonlinear_state - obj.g(obj.particles{j}.prev_nonlinear_state);
%                 cov_qy = obj.Cov_qxy(obj.n_linear+1:end,obj.n_linear+1:end);
%                 Gx = obj.Gx(obj.particles{j}.nonlinear_state);                
%                 obj.particles{j} = obj.particles{j}.linear_update(zeros(size(Gx,1), 1), Gx, cov_qy, z_dyn);                      

                
                
            end            
        end    
        
        function obj = filter2(obj, z)
            weights_sum = 0;
            
            for (j = 1:length(obj.particles))
                % 1. Non-linear state Particle filter prediction
                % Draw samples from the non-linear proposal density/distribution
                % Note that when drawing samples from the non-linear
                % distribution the uncertainty of the linear state has to be
                % considered by drawing a sample of the linear state and use
                % this when drawing a sample from the conditional distribution
                % for the non-linear propgation model            
                y_prev = obj.particles{j}.nonlinear_state;
                % Compute resulting noise covariance of particle by
                % inflating the noise with the uncertainty in the linear state                
                cov_nonlinear = obj.Gx(y_prev) * obj.particles{j}.P * obj.Gx(y_prev)' + obj.Cov_qxy(obj.n_linear+1:end,obj.n_linear+1:end);                                
                % Draw noise sample
                qy = chol(cov_nonlinear, 'lower') * randn(obj.n_nonlinear,1);            
                % Compute new non-linear state
                y = obj.g(y_prev) + obj.Gx(y_prev) * obj.particles{j}.x + qy;
                obj.particles{j}.nonlinear_state = y;
                obj.particles{j}.prev_nonlinear_state = y_prev;            
            
                % 2. Linear state Kalman filter update with change in the non-linear state
                % The sample that was just drawn for the non-linear state is
                % used in an update step of the Kalman filter for the linear
                % state. This thus requires that an dynamics-based measurement
                % model is constructed for the linear state, relating a
                % previous non-linear state to the current new non-linear state
                % Note that this can be seen as using the motion model of the
                % non-linear state as a measurement model for the linear state.
                % y[k] = g(y[k-1]) + Gx(y[k-1]) * x[k-1] + qy[k-1]
                % y[k] - g(y[k-1]) = Gx(y[k-1]) * x[k-1] + qy[k-1]            
                z_dyn = obj.particles{j}.nonlinear_state - obj.g(obj.particles{j}.prev_nonlinear_state);
                cov_qy = obj.Cov_qxy(obj.n_linear+1:end,obj.n_linear+1:end);
                Gx = obj.Gx(obj.particles{j}.nonlinear_state);                
                obj.particles{j} = obj.particles{j}.linear_update(zeros(size(Gx,1), 1), Gx, cov_qy, z_dyn);            
            
                % 3. Linear state Kalman filter prediction
                % Apply the linear motion model to predict/propagate the linear
                % state and covariance                            
                cov_qx = obj.Cov_qxy(1:obj.n_linear,1:obj.n_linear);
                f = obj.f(obj.particles{j}.nonlinear_state);
                Fx = obj.Fx(obj.particles{j}.nonlinear_state);
                obj.particles{j} = obj.particles{j}.linear_predict(f, Fx, cov_qx);                
            
                % 4. Non-linear state Particle filter update
                % Use the measurement model likelihood function to update the
                % weight of each particle. Note that if the measurement model
                % depends on the linear state(s), these will have to be
                % marginalized out when calculating this likelihood.
                y = obj.particles{j}.nonlinear_state;
                likelihood = obj.measurement_nonlinear_pdf(z, y);                                                
                obj.particles{j}.weight = obj.particles{j}.weight * likelihood;
                weights_sum = weights_sum + obj.particles{j}.weight;
                
                % 5. Linear state Kalman filter measurement update
                % Perform a measurement update of the linear states using the
                % measurement model where the non-linear states are given by
                % the particle state (viewed as a known constant/deterministic).                
                h = obj.h(obj.particles{j}.nonlinear_state);
                Hx = obj.Hx(obj.particles{j}.nonlinear_state);     
                if (any(Hx))
                    obj.particles{j} = obj.particles{j}.linear_update(h, Hx, obj.Cov_r, z);                            
                end
            end
            
            % Normalize weights
            for (j = 1:length(obj.particles))
                obj.particles{j}.weight = obj.particles{j}.weight / weights_sum;
            end                        
                
            % Check for resampling need            
            N = length(obj.particles);
            if (obj.getEffectiveNumberOfParticles() < N/4)
                disp('Resampling');
                obj = obj.resample();
            end
        end        
        
        function nonlinear_states = getNonlinearStates(obj)
            nonlinear_states = zeros(obj.n_nonlinear, length(obj.particles));
            for (j = 1:length(obj.particles))
                nonlinear_states(:,j) = obj.particles{j}.nonlinear_state;  
            end
        end
        
        function linear_states = getLinearStates(obj)
            linear_states = zeros(obj.n_nonlinear, length(obj.particles));
            for (j = 1:length(obj.particles))
                linear_states(:,j) = obj.particles{j}.x;  
            end
        end        
        
        function linear_gaussians = getLinearGaussians(obj)
            linear_gaussians = {};
            for (j = 1:length(obj.particles))
                linear_gaussians{j} = Gaussian(obj.particles{j}.x, obj.particles{j}.P);  
            end
        end        
                
        
        function combined_states = getCombinedStates(obj)
            combined_states = zeros(obj.n_linear + obj.n_nonlinear, length(obj.particles));
            for (j = 1:length(obj.particles))                
                combined_states(:,j) = [obj.particles{j}.x;  
                                   obj.particles{j}.nonlinear_state];
            end
        end
        
        function [states, weights] = getStatesAndWeights(obj)
            states = zeros(obj.n_linear + obj.n_nonlinear, length(obj.particles));
            weights = zeros(1, length(obj.particles));
            for (j = 1:length(obj.particles))                
                states(:,j) = [obj.particles{j}.x;  
                                   obj.particles{j}.nonlinear_state];
                weights(:,j) = obj.particles{j}.weight;
            end
        end        
    end        
    
    methods %(Access = private)
        function Neff = getEffectiveNumberOfParticles(obj)
            % Use this to identify when to resample (when this number goes below N/4)
            sum_of_squared_weights = 0;
            for (j = 1:length(obj.particles))
                sum_of_squared_weights = sum_of_squared_weights + obj.particles{j}.weight^2;
            end
            Neff = 1 / sum_of_squared_weights;                
        end
        
        function obj = resample(obj)
            % Draw N samples with replacement from the current set of
            % particles with the probability of drawing a particular
            % particle defined according to its' weight.
            % We draw the samples by drawing a number from a uniform
            % distribution between 0 to 1 and then using a cumulated sum
            % lookup table to find the index of the according particle.
            new_particles = {};
                        
            %resampling_wheel = cumsum(obj.weights);
            resampling_wheel = zeros(1,length(obj.particles));
            resampling_wheel(1) = obj.particles{1}.weight;
            for (j = 2:length(obj.particles))                
                resampling_wheel(j) = resampling_wheel(j-1) + obj.particles{j}.weight;
            end
                
            uniform_samples = rand(1, length(obj.particles));
            for (i = 1:length(obj.particles))                 
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
                new_particles{i} = obj.particles{old_particle_idx};
            end
            obj.particles = new_particles;
            
            % Set all the weights and uniformly
            weight = 1 / length(obj.particles);
            for (j = 1:length(obj.particles))                
                obj.particles{j}.weight = weight;
            end            
        end     
        
        function obj = low_variance_resample(obj)
            % Draw N samples from the current set as a sequential stochastic
            % process with the probability of drawing a particular particle
            % defined according to its' weight. Note that this is different
            % from the normal resampling strategy where samples are drawn
            % independently.
            %
            % The resampling step is a probabilistic implementation of the
            % Darwinian idea of survival of the fittest: It refocuses the
            % particle set to regions in state space with high posterior
            % probability
            %
            % We draw the samples by drawing ONE random number, r, from a
            % uniform distribution between 0 to 1/N and then using a cumulated sum
            % lookup table to find the index of the according particle.
            new_particles = {};
                        
            %resampling_wheel = cumsum(obj.weights);
            resampling_wheel = zeros(1,length(obj.particles));
            resampling_wheel(1) = obj.particles{1}.weight;
            for (j = 2:length(obj.particles))                
                resampling_wheel(j) = resampling_wheel(j-1) + obj.particles{j}.weight;
            end
            
            N = length(obj.particles);
            r = rand(1, 1) / N;  
            u = r;
            % OBS. That the below can be optimized (made more efficiently)
            % according to the algorithm "Algorithm Low variance sampler"
            % in Probabilistic Robotics by Sebastian Thrun
            for (i = 1:N)                
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
                new_particles{i} = obj.particles{old_particle_idx};
                u = u + 1/N;
            end
            obj.particles = new_particles;
            
            % Set all the weights and uniformly
            weight = 1 / length(obj.particles);
            for (j = 1:length(obj.particles))                
                obj.particles{j}.weight = weight;
            end  
        end        
        
        function obj = resample_other_strategies(obj)
%            switch resampling_strategy
%                case 'multinomial_resampling'
%                   with_replacement = true;
%                   idx = randsample(1:Ns, Ns, with_replacement, wk);
%             %{
%                   THIS IS EQUIVALENT TO:
%                   edges = min([0 cumsum(wk)'],1); % protect against accumulated round-off
%                   edges(end) = 1;                 % get the upper edge exact
%                   % this works like the inverse of the empirical distribution and returns
%                   % the interval where the sample is to be found
%                   [~, idx] = histc(sort(rand(Ns,1)), edges);
%             %}
%                case 'systematic_resampling'
%                   % this is performing latin hypercube sampling on wk
%                   edges = min([0 cumsum(wk)'],1); % protect against accumulated round-off
%                   edges(end) = 1;                 % get the upper edge exact
%                   u1 = rand/Ns;
%                   % this works like the inverse of the empirical distribution and returns
%                   % the interval where the sample is to be found
%                   [~, idx] = histc(u1:1/Ns:1, edges);
%                % case 'regularized_pf'      TO BE IMPLEMENTED
%                % case 'stratified_sampling' TO BE IMPLEMENTED
%                % case 'residual_sampling'   TO BE IMPLEMENTED
%                otherwise
%                   error('Resampling strategy not implemented')
%             end;
        end
    end
end