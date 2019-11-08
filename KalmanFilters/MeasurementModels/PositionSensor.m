function [h, Hx, Hr] = PositionSensor(idx_x, idx_y)

    % Global position sense    
    % See "5.6.1 Measurement models" from Course ChM015x
    
    % State vector assumed to include
    %   x = [ x, y, ... ]
    % Measurement vector:
    %   z = [ zx, zy ]
    % Noise vector:
    %   r = [ r_x, r_ry ]    
    
    % Measurement model
    %   z_x = x + r_x    
    %   z_y = y + r_y
    
    % x is state    
    % r is noise vector
    h = @(x, r) [ ...
        x(idx_x) + r(1)
        x(idx_y) + r(2)
    ];

    % Jacobians of motion model
    Hx = @(x, r) [ ...
        zeros(1, idx_x-1), 1, zeros(1, size(x,1)-idx_x)
        zeros(1, idx_y-1), 1, zeros(1, size(x,1)-idx_y)
    ];

    Hr = @(x, r) [ ...
        1, 0
        0, 1
    ];

end