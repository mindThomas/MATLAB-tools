function [f, Fx, Fu] = CoordinatedTurnModel_Discrete2(ts)
    % Coordinated Turn model from Probabilistic Robotics
    % See Equation 5.9 in Probabilistic Robotics by Sebastic Thrun
    % This model is slightly different than the model presented in Course
    % ChM015x
    %
    % Exact solution to discrete time kinematic coordinated turn model    
        
    % State vector:
    %   x = [ x, y, phi ]
    % Input vector:
    %   u = [ v, omega ] 
    
    % Difference equations:
    %   x[k] = x[k-1] - v[k-1]/omega[k-1] * sin(phi[k-1]) + v[k-1]/omega[k-1] * sin(phi[k-1] + omega[k-1]*ts)
    %   y[k] = y[k-1] + v[k-1]/omega[k-1] * cos(phi[k-1]) - v[k-1]/omega[k-1] * cos(phi[k-1] + omega[k-1]*ts) 
    %   phi[k] = phi[k-1] + ts*omega[k-1]    
       
    % x is the previous state, x[k-1]
    % u is the previous inputs, u[k-1]
    % q is the discrete zero-order hold noise vector, q[k-1]
    f = @(x, u) x + [ ...
            -u(1)/u(2)*sin(x(3)) + u(1)/u(2)*sin(x(3)+ts*u(2))
            u(1)/u(2)*cos(x(3)) - u(1)/u(2)*cos(x(3)+ts*u(2))
            ts*u(2)
        ];

    % Jacobians of motion model
    Fx = @(x, u) eye(3) + [ ...
        0, 0, u(1)/u(2) * (-cos(x(3)) + cos(x(3)+ts*u(2)))
        0, 0, u(1)/u(2) * (-sin(x(3)) + sin(x(3)+ts*u(2)))
        0, 0, 0        
    ];

    Fu = @(x, u) [ ...
        -1/u(2)*sin(x(3)) + 1/u(2)*sin(x(3)+ts*u(2)), u(1)/u(2)^2*sin(x(3)) - u(1)/u(2)^2*sin(x(3)+ts*u(2)) + u(1)/u(2)*ts*cos(x(3)+ts*u(2))
        1/u(2)*cos(x(3)) - 1/u(2)*cos(x(3)+ts*u(2)), -u(1)/u(2)^2*cos(x(3)) + u(1)/u(2)^2*cos(x(3)+ts*u(2)) + u(1)/u(2)*ts*sin(x(3)+ts*u(2))
        0, ts
    ];
    
end