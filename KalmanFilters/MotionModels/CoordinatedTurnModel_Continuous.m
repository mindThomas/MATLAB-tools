function [f, Fx, Fu, Fq] = CoordinatedTurnModel_Continuous()

    % Continuous time coordinated turn model
    % A kinematic coordinated turn model
    %   Velocity, v, is a Wiener process, that is acceleration is a Gaussian
    %     noise random variable
    %   Heading, phi, is described by a continuous velocity model
    
    % See "5.5.1 A coordinated turn model" from Course ChM015x
    
    % State vector:
    %   x = [ x, y, v, phi, omega ]
    % Input vector:
    %   u = []
    % Noise vector:
    %   q = [ q_v, q_omega ]    
    
    % Differential equations:
    %   dx = v * cos(phi)
    %   dy = v * sin(phi)
    %   dv = 0 + q_v
    %   dphi = omega
    %   domega = 0 + q_omega
    
    % x is state
    % u is inputs
    % q is noise vector
    f = @(x, u, q) [ ...
        x(3) * cos(x(4))
        x(3) * sin(x(4))
        q(1)
        x(5)
        q(2)
    ];

    % df(1)/dx(3) = cos(x(4))
    % df(1)/dx(4) = x(3)*dcos(g)/dg * 1 = -x(3)*sin(x(4))
    % df(2)/dx(3) = sin(x(4))
    % df(2)/dx(4) = x(3)*cos(x(4))
    
    % Jacobians of motion model
    Fx = @(x, u, q) [ ...
        0, 0, cos(x(4)), -x(3)*sin(x(4)), 0
        0, 0, sin(x(4)), x(3)*cos(x(4)), 0
        0, 0, 0, 0, 0
        0, 0, 0, 0, 1
        0, 0, 0, 0, 0
    ];

    Fu = @(x, u, q) [ ...
    ];

    Fq = @(x, u, q) [ ...
        0, 0
        0, 0
        1, 0
        0, 0
        0, 1
    ];

end