function [f, Fx, Fu, Fq] = CoordinatedTurnModel_Discrete(ts)

    % Exact solution to discrete time kinematic coordinated turn model    
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
    
    % Difference equations:
    %   x[k] = x[k-1] + 2*v[k-1]/omega[k-1] * sin(omega[k-1]*ts/2) * cos(phi[k-1] + omega[k-1]*ts/2)
    %   y[k] = y[k-1] + 2*v[k-1]/omega[k-1] * sin(omega[k-1]*ts/2) * sin(phi[k-1] + omega[k-1]*ts/2)
    %   v[k] = v[k-1] + q_v[-1]
    %   phi[k] = phi[k-1] + ts*omega[k-1]
    %   omega[k] = omega[k-1] + q_omega[k-1]    
    
    % Smart reciprocal
    %rec = @(x) prod(1/x*(abs(x)>eps), 'omitnan')*(abs(x)>eps);
    rec = @(x) 1/x;
    
    % x is the previous state, x[k-1]
    % u is the previous inputs, u[k-1]
    % q is the discrete zero-order hold noise vector, q[k-1]
    f = @(x, u, q) x + [ ...
         	2*x(3)*rec(x(5)) * sin(x(5)*ts/2) * cos(x(4) + x(5)*ts/2) + ts*x(3)*cos(x(4))*(abs(x(5))<=eps)
            2*x(3)*rec(x(5)) * sin(x(5)*ts/2) * sin(x(4) + x(5)*ts/2) + ts*x(3)*sin(x(4))*(abs(x(5))<=eps)
            q(1)
            ts*x(5)
            q(2)
        ];

    % dx(1)/dx(1) = 0
    % dx(1)/dx(3) = 2/x(5) * sin(x(5)*ts/2) * cos(x(4) + x(5)*ts/2) + 
    %               2*x(3)/x(5) * 0 * cos(x(4) + x(5)*ts/2) + 
    %               2*x(3)/x(5) * sin(x(5)*ts/2) * 0
    % dx(1)/dx(3) = 2/x(5) * sin(x(5)*ts/2) * cos(x(4) + x(5)*ts/2)
    
    % dx(1)/dx(4) = 0 * sin(x(5)*ts/2) * cos(x(4) + x(5)*ts/2) + 
    %               2*x(3)/x(5) * 0 * cos(x(4) + x(5)*ts/2) +
    %               2*x(3)/x(5) * sin(x(5)*ts/2) * ...
    % dcos(g)/dg * d(x(4) + x(5)*ts/2)/dx(4)
    % = -sin(x(4) + x(5)*ts/2)
    % dx(1)/dx(4) = -2*x(3)/x(5) * sin(x(5)*ts/2) * sin(x(4) + x(5)*ts/2)
    
    % dx(1)/dx(5) = -2*x(3)/x(5)^2 * sin(x(5)*ts/2) * cos(x(4) + x(5)*ts/2) + 
    %               2*x(3)/x(5) * ts/2*cos(x(5)*ts/2) * cos(x(4) + x(5)*ts/2) +
    %               2*x(3)/x(5) * sin(x(5)*ts/2) * ts/2*cos(x(4) + x(5)*ts/2)    
    % dx(1)/dx(5) = -2*x(3)/x(5)^2 * sin(x(5)*ts/2) * cos(x(4) + x(5)*ts/2) + 
    %               x(3)/x(5) * ts*cos(x(5)*ts/2) * cos(x(4) + x(5)*ts/2) +
    %               -x(3)/x(5) * sin(x(5)*ts/2) * ts*sin(x(4) + x(5)*ts/2)        
    
    % Jacobians of motion model
    Fx = @(x, u, q) eye(5) + [ ...
        0, 0, 2*rec(x(5))*sin(x(5)*ts/2)*cos(x(4)+x(5)*ts/2) + ts*cos(x(4))*(abs(x(5))<=eps),   -2*x(3)*rec(x(5))*sin(x(5)*ts/2)*sin(x(4)+x(5)*ts/2) - ts*x(3)*sin(x(4))*(abs(x(5))<=eps),         x(3)*rec(x(5))^2 * (x(5)*ts*cos(x(5)*ts/2)*cos(x(5)*ts/2 + x(4)) - x(5)*ts*sin(x(5)*ts/2)*sin(x(5)*ts/2+x(4)) - 2*sin(x(5)*ts/2)*cos(x(5)*ts/2 + x(4)) )
        0, 0, 2*rec(x(5))*sin(x(5)*ts/2)*sin(x(4)+x(5)*ts/2) + ts*sin(x(4))*(abs(x(5))<=eps),    2*x(3)*rec(x(5)) * sin(x(5)*ts/2) * cos(x(4) + x(5)*ts/2) + ts*x(3)*cos(x(4))*(abs(x(5))<=eps),   x(3)*rec(x(5))^2 * (x(5)*ts*sin(x(5)*ts/2)*cos(x(5)*ts/2 + x(4)) + x(5)*ts*cos(x(5)*ts/2)*sin(x(5)*ts/2+x(4)) - 2*sin(x(5)*ts/2)*sin(x(5)*ts/2 + x(4)) )
        0, 0, 0, 0, 0
        0, 0, 0, 0, ts
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