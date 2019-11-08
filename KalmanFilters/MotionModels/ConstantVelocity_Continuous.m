function [f, Fx, Fu, Fq] = ConstantVelocity_Continuous()

    % Continuous time kinematic constant velocity model    
    %   Velocity, v, is a Wiener process, that is acceleration is a Gaussian
    %     noise random variable    
    
    % See "5.3.1 Discretizing linear models: the transition matrix" from Course ChM015x
    
    % State vector:
    %   x = [ p, v ]
    % Input vector:
    %   u = []
    % Noise vector:
    %   q = [ q_v ]
    
    % Differential equations:
    %   dp = v    
    %   dv = 0 + q_v    
    
    % x is state
    % u is inputs
    % q is noise vector
    f = @(x, u, q) [ ...
        x(2)        
        q(1)        
    ];

    % Jacobians of motion model
    Fx = @(x, u, q) [ ...
        0, 1
        0, 0
    ];

    Fu = @(x, u, q) [ ...
    ];

    Fq = @(x, u, q) [ ...
        0
        1
    ];

end