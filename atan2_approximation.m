function res = atan2_approximation(y, x)
    % http://pubs.opengroup.org/onlinepubs/009695399/functions/atan2.html
    % Volkan SALMA : https://gist.github.com/volkansalma/2972237   
    abs_y = abs(y);

    if (x < 0)    
        r = (x + abs_y) ./ (abs_y - x);
        angle = 3 * pi / 4; % three_quarter_pi;    
    else    
        r = (x - abs_y) ./ (x + abs_y);
        angle = pi / 4; % one_quarter_pi;
    end

    angle = angle + (0.1963 .* r .* r - 0.9817) .* r;

    if (y < 0)    
        res = -angle; % negate if in quad III or IV    
        return;
    end
        
    res = angle;
end