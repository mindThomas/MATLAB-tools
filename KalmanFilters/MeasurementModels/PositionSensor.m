function [h, Hx, Hr] = PositionSensor(idx_p)

    % Global position sense    
    % See "5.6.1 Measurement models" from Course ChM015x
    
    % State vector assumed to include
    %   x = [ p, ... ]
    % Measurement vector:
    %   z = [ zp ]
    % Noise vector:
    %   r = [ r_p ]    
    
    % Measurement model
    %   z_p = p + r_p    
    
    % x is state    
    % r is noise vector
    h = @(x, r) [ ...
        x(idx_p) + r(1)        
    ];

    % Jacobians of motion model
    Hx = @(x, r) [ ...
        zeros(1, idx_p-1), 1, zeros(1, size(x,1)-idx_p)        
    ];

    Hr = @(x, r) [ ...
        1        
    ];

end