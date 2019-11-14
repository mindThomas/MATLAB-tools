function [h, Hx, Hr] = RangeBearingSensor(idx_x, idx_y)

    % Global position sense    
    % See "5.6.1 Measurement models" from Course ChM015x
    
    % State vector assumed to include
    %   x = [ x, y, ... ]
    % Measurement vector:
    %   z = [ z_r, z_phi ]
    % Noise vector:
    %   r = [ r_r, r_phi ]    
    
    % Measurement model
    %   z_r = sqrt( x^2 + y^2 ) + r_r    
    %   z_phi = atan(y / x) + r_phi
    
    % x is state    
    % r is noise vector
    h = @(x, r) [ ...
        sqrt( x(idx_x)^2 + x(idx_x)^2 ) + r(1)
        atan2(x(idx_y), x(idx_x)) + r(2)
    ];

    % Jacobians of motion model
    Hx = @(x, r) ...
        [zeros(1, idx_x-1), x(idx_x) / sqrt( x(idx_x)^2 + x(idx_y)^2), zeros(1, size(x,1)-idx_x)
         zeros(1, idx_x-1), -x(idx_y) / sqrt( x(idx_x)^2 + x(idx_y)^2), zeros(1, size(x,1)-idx_x)] + ...
        [zeros(1, idx_y-1), x(idx_y) / sqrt( x(idx_x)^2 + x(idx_y)^2), zeros(1, size(x,1)-idx_y)
         zeros(1, idx_y-1), x(idx_x) / sqrt( x(idx_x)^2 + x(idx_y)^2), zeros(1, size(x,1)-idx_y)];

    Hr = @(x, r) [ ...
        1, 0
        0, 1
    ];

end