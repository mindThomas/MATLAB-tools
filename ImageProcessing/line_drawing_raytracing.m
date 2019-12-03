function line_drawing_raytracing(width, height, center_x, center_y, angle, length)
    % Bresenham's line tracing algorithm
    % https://www.coursera.org/lecture/motion-planning-self-driving-cars/lesson-2-populating-occupancy-grids-from-lidar-scan-data-part-2-VcH67
    % From https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
    im = zeros(height, width);
    
    x0 = round(center_x);
    y0 = round(center_y);    
    
    x1 = center_x + length * cos(angle);
    y1 = center_y + length * sin(angle);
    x1 = round(x1);
    y1 = round(y1);

    deltax = x1 - x0;
    deltay = y1 - y0;
    if (abs(deltax) > abs(deltay))
        deltaerr = abs(deltay / deltax); % Assume deltax != 0 (line is not vertical),
        % note that this division needs to be done in a way that preserves the fractional part
        error = 0.0; % No error at start
        y = y0;
        for (x = x0:x1)
            im(y+1,x+1) = 1;
            error = error + deltaerr;
            if (error >= 0.5)
                y = y + sign(deltay) * 1;
                error = error - 1.0;
            end
        end
    else
        deltaerr = abs(deltax / deltay); % Assume deltax != 0 (line is not vertical),
        % note that this division needs to be done in a way that preserves the fractional part
        error = 0.0; % No error at start
        x = x0;
        for (y = y0:y1)
            im(y,x) = 1;
            error = error + deltaerr;
            if (error >= 0.5)
                x = x + sign(deltax) * 1;
                error = error - 1.0;
            end
        end
    end

    imshow(im(end:-1:1,:));
end