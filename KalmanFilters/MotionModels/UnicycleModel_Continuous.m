function [f, Fx, Fu, Fq] = UnicycleModel_Continuous(u)

    % Continuous time kinematic unicycle model    
    % See also the Coordinated turn model   
    
    % State vector:
    %   x = [ x, y, psi ]
    % Input vector:
    %   u = [ v, omega ]
    % Noise vector:
    %   q = [ q_v, q_omega ]    
    
    % Differential equations:
    %   dx = (v + q_v) * cos(psi)
    %   dy = (v + q_v) * sin(psi)    
    %   dpsi = omega + q_omega    
    
    % x is state
    % u is inputs
    % q is noise vector
    f = @(x, u, q) [ ...
        (u(1) + q(1)) * cos(x(3))
        (u(1) + q(1)) * sin(x(3))
        x(2) + q(2)
    ];
    
    % Jacobians of motion model
    Fx = @(x, u, q) [ ...
        0, 0, -(u(1) + q(1)) * sin(x(3))
        0, 0, (u(1) + q(1)) * cos(x(3))
        0, 1, 0
    ];

    Fu = @(x, u, q) [ ...
        cos(x(3)), 0
        sin(x(3)), 0
        0, 1
    ];

    Fq = @(x, u, q) [ ...
        cos(x(3)), 0
        sin(x(3)), 0
        0, 1
    ];

end