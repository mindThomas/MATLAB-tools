function [f, Fx, Fu, Fq] = ConstantAcceleration2D_Continuous()

    % Continuous time kinematic constant acceleration model    
    %   Acceleration, a, is a Wiener process, that is the derivative of
    %   the acceleration is a Gaussian noise random variable    
    
    % See "5.3.1 Discretizing linear models: the transition matrix" from Course ChM015x
    
    % State vector:
    %   x = [ x, y, vx, vy, ax, ay ]
    % Input vector:
    %   u = []
    % Noise vector:
    %   q = [ q_ax, q_ay ]    
    
    % Differential equations:
    %   dx = vx
    %   dy = vy
    %   dvx = ax
    %   dvy = ay
    %   dax = 0 + q_ax
    %   day = 0 + q_ay
    
    % x is state
    % u is inputs
    % q is noise vector
    f = @(x, u, q) [ ...
        x(3)
        x(4)
        x(5)
        x(6)
        q(1)
        q(2)
    ];

    % Jacobians of motion model
    Fx = @(x, u, q) [ ...
        0, 0, 1, 0, 0, 0
        0, 0, 0, 1, 0, 0
        0, 0, 0, 0, 1, 0
        0, 0, 0, 0, 0, 1
        0, 0, 0, 0, 0, 0
        0, 0, 0, 0, 0, 0
    ];

    Fu = @(x, u, q) [ ...
    ];

    Fq = @(x, u, q) [ ...
        0, 0
        0, 0
        0, 0
        0, 0
        1, 0
        0, 1
    ];

end