function [f, Fx, Fu, Fq] = ConstantAcceleration_Continuous()

    % Continuous time kinematic constant acceleration model    
    %   Acceleration, a, is a Wiener process, that is the derivative of
    %   the acceleration is a Gaussian noise random variable    
    
    % See "5.3.1 Discretizing linear models: the transition matrix" from Course ChM015x
    
    % State vector:
    %   x = [ p, v, a ]
    % Input vector:
    %   u = []
    % Noise vector:
    %   q = [ q_a ]    
    
    % Differential equations:
    %   dp = v    
    %   dv = a
    %   da = 0 + q_a
    
    % x is state
    % u is inputs
    % q is noise vector
    f = @(x, u, q) [ ...
        x(2)
        x(3)
        q(1)        
    ];

    % Jacobians of motion model
    Fx = @(x, u, q) [ ...
        0, 1, 0
        0, 0, 1
        0, 0, 0
    ];

    Fu = @(x, u, q) [ ...
    ];

    Fq = @(x, u, q) [ ...
        0
        0
        1
    ];

end