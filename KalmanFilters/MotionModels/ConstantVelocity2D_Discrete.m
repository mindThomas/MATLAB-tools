function [f, Fx, Fu, Fq] = ConstantVelocity2D_Discrete(ts)

    % Dicrete time kinematic constant velocity model    
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
    %   x[k] = x[k-1] + ts * vx[k-1]
    %   y[k] = y[k-1] + ts * vy[k-1]
    %   vx[k] = vx[k-1] + q_vx[k-1]
    %   vy[k] = vy[k-1] + q_vy[k-1]
    
    % x is state
    % u is inputs
    % q is noise vector
    f = @(x, u, q) x + [ ...
        ts * x(3)
        ts * x(4)
        q(1)
        q(2)
    ];

    % Jacobians of motion model
    Fx = @(x, u, q) eye(4) + [ ...
        0, 0, ts, 0
        0, 0, 0, ts
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