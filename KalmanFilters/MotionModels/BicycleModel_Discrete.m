function [f, Fx, Fu, Fq] = BicycleModel_Discrete(L, dt)

    % Discrete time kinematic bicycle model with origin in the rear axis
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
    
    % Difference equations:
    % dpsi = v * 1/L * tan(rho)
    %   x[k] = x[k-1] + v / dpsi * (sin(psi + dt*psidot) - sin(psi))
    %   y[k] = y[k-1] v / dpsi * (-cos(psi + dt*psidot) + cos(psi))
    %   The above is similar to Coodinated Turn, but still fairly different
    %   psi[k] = psi[k-1] + v * 1/L * tan(rho)
    %   v[k] = v[k-1] + q_v[k]
    %   rho[k] = rho[k-1] + q_rho[k]
    
    % x is state
    % u is inputs
    % q is noise vector
    f = @(x, u, q) [ ...
    ];

    % Jacobians of motion model
    Fx = @(x, u, q) [ ...
    ];

    Fu = @(x, u, q) [ ...
    ];

    Fq = @(x, u, q) [ ...
    ];

end