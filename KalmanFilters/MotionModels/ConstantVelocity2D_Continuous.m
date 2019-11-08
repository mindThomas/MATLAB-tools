function [f, Fx, Fu, Fq] = ConstantVelocity2D_Continuous()

    % Continuous time kinematic constant velocity model    
    %   Velocity, v, is a Wiener process, that is acceleration is a Gaussian
    %     noise random variable    
    
    % See "5.3.1 Discretizing linear models: the transition matrix" from Course ChM015x
    
    % State vector:
    %   x = [ x, y, vx, vy ]
    % Input vector:
    %   u = []
    % Noise vector:
    %   q = [ q_vx, q_vy ]    
    
    % Differential equations:
    %   dx = vx
    %   dy = vy
    %   dvx = 0 + q_vx
    %   dvy = 0 + q_vy    
    
    % x is state
    % u is inputs
    % q is noise vector
    f = @(x, u, q) [ ...
        x(3)
        x(4)
        q(1)        
        q(2)
    ];

    % Jacobians of motion model
    Fx = @(x, u, q) [ ...
        0, 0, 1, 0
        0, 0, 0, 1
        0, 0, 0, 0
        0, 0, 0, 0        
    ];

    Fu = @(x, u, q) [ ...
    ];

    Fq = @(x, u, q) [ ...
        0, 0
        0, 0
        1, 0
        0, 1
    ];

end