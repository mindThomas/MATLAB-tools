classdef BinaryBayesFilter
    % Implementation of a Binary Bayes filter
    % See Table 4.2 in "Probabilistic Robotics" by Sebastian Thrun
    properties (Access = private)
        %p % p(x) = p(x=1) : probability of the state being 1
        l % log odds of the state = ratio of the probability divided by the probability of its negate
        l0 % initial log odds
        h_inv % inverse measurement distribution, p(x[k] | z[k]), on the form @(z, x) ...
              % where the state is binary and can thus only take the values
              % x[k]=0 and x[k]=1
    end
    methods
        function obj = BinaryBayesFilter(h_inv, init_probability)           
            obj.h_inv = h_inv;
            p0 = init_probability;
            obj.l0 = logOdds(p0);
        end             
        
        % The Bayes filter for updating beliefs
        % in log odds representation is computationally elegant. It avoids 
        % truncation problems that arise for probabilities close to 0 or 1.
        function update(obj, z)
            % The binary Bayes filter uses an inverse measurement model
            % p(x | zt), instead of the familiar forward model p(zt | x).
            % The inverse measurement model specifies a distribution over
            % the (binary) state variable as a function of the measurement z.
            % See equation 4.15 to 4.20 in "Probabilistic Robotics"
            p_x_given_z = h_inv(z, 1); % probability of x=1 given the measurement            
            obj.l = obj.l + logOdds(p_x_given_z) - obj.l0;
        end
        
        function odds = odds(p)
            % The odds of a state x is defined as the ratio of the
            % probability of this event divided by the probability of its negate
            odds = p / (1 - p);
        end   
        
        function l = logOdds(p)
            % For binary states the belief (posterior of the state) is
            % commonly implemented as a log odds ratio.
            % Log odds assume values from -inf to inf
            l = log( p / (1 - p) );
        end 
        
        function p = getProbability(obj)
            % Convert log odds to probability probability of the state being 1
            % p(x) = p(x=1)
            p = 1 - 1 / (1+exp(obj.l));
        end
    end
end