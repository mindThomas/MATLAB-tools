function [f, Fx, Fu, Fq] = ConstantVelocity_Discrete(ts)

    % Dicrete time kinematic constant velocity model    
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
    %   p[k] = p[k-1] + ts * v[k-1]
    %   v[k] = v[k-1] + q_v[k-1]
    
    % x is state
    % u is inputs
    % q is noise vector
    f = @(x, u, q) x + [ ...
        ts * x(2)        
        q(1)        
    ];

    % Jacobians of motion model
    Fx = @(x, u, q) eye(2) + [ ...
        0, ts
        0, 0
    ];

    Fu = @(x, u, q) [ ...
    ];

    Fq = @(x, u, q) [ ...
        0
        1
    ];

end