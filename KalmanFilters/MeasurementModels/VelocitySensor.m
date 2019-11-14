function [h, Hx, Hr] = VelocitySensor(idx_vx, idx_vy)

    % Global position sense    
    % See "5.6.1 Measurement models" from Course ChM015x
    
    % State vector assumed to include
    %   x = [ x, y, vx, vy, ... ]
    % Measurement vector:
    %   z = [ z_v ]
    % Noise vector:
    %   r = [ r_v ]    
    
    % Measurement model
    %   z_v = sqrt( vx^2 + vy^2 )
    
    % x is state    
    % r is noise vector
    h = @(x, r) [ ...
        sqrt( x(idx_vx)^2 + x(idx_vy)^2) + r(1)        
    ];

    % Jacobians of motion model
    Hx = @(x, r) ...
        [zeros(1, idx_vx-1), x(idx_vx) / sqrt( x(idx_vx)^2 + x(idx_vy)^2), zeros(1, size(x,1)-idx_vx)] + ...
        [zeros(1, idx_vy-1), x(idx_vy) / sqrt( x(idx_vx)^2 + x(idx_vy)^2), zeros(1, size(x,1)-idx_vy)];    

    Hr = @(x, r) [ ...
        1        
    ];

end