function [f, Fx, Fu, Fq] = BicycleModel_Continuous(L)

    % Continuous time kinematic bicycle model with origin in the rear axis
    %   Velocity, v, and steering angle, rho, are Wiener processes,
    %     that is acceleration is a Gaussian noise random variable   
    % Constants:
    %   L : the distance from front to rear wheels
    % See:
    %   - https://www.ri.cmu.edu/pub_files/2009/2/Automatic_Steering_Methods_for_Autonomous_Automobile_Path_Tracking.pdf
    %   - https://www.coursera.org/lecture/intro-self-driving-cars/lesson-2-the-kinematic-bicycle-model-Bi8yE
    
    % State vector:
    %   x = [ x, y, psi, v, rho ]
    % Input vector:
    %   u = []
    % Noise vector:
    %   q = [ q_v, q_rho ]    
    
    % Differential equations:
    %   dx = v * cos(psi)
    %   dy = v * sin(psi)   
    %   dpsi = v * 1/L * tan(rho)
    %   dv = 0 + q_v
    %   drho = 0 + q_rho
    
    % x is state
    % u is inputs
    % q is noise vector
    f = @(x, u, q) [ ...
        x(4) * cos(x(3))
        x(4) * sin(x(3))
        x(4) * 1/L * tan(x(5))
        q(1)
        q(2)
    ];

    % Jacobians of motion model
    Fx = @(x, u, q) [ ...
        0, 0, -x(4)*sin(x(3)), cos(x(3)), 0
        0, 0, x(4)*cos(x(3)), sin(x(3)), 0
        0, 0, 0, 1/L*tan(x(5)), x(4)*1/L*cos(x(5))^2
        0, 0, 0, 0, 0     
        0, 0, 0, 0, 0
    ];

    Fu = @(x, u, q) [ ...
    ];

    Fq = @(x, u, q) [ ...
        0, 0
        0, 0
        0, 0
        1, 0
        0, 1
    ];

end