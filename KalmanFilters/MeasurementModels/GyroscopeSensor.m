function [h, Hx, Hr] = GyroscopeSensor(idx_omega)

    % Global position sense    
    % See "5.6.1 Measurement models" from Course ChM015x
    
    % State vector assumed to include
    %   x = [ x, y, vx, vy, phi, omega, ... ]
    % Measurement vector:
    %   z = [ z_omega ]
    % Noise vector:
    %   r = [ r_omega ]    
    
    % Measurement model
    %   z_omega = omega + r_omega
    
    % x is state    
    % r is noise vector
    h = @(x, r) [ ...
        x(idx_omega) + r(1)        
    ];

    % Jacobians of motion model
    Hx = @(x, r) [ ...
        zeros(1, idx_omega-1), 1, zeros(1, size(x,1)-idx_omega)        
    ];

    Hr = @(x, r) [ ...
        1        
    ];

end