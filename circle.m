function h = circle(x, y, r, varargin)        
    draw_resolution = deg2rad(1);
    theta = 0:draw_resolution:2*pi;
    xx = r * cos(theta) + x;
    yy = r * sin(theta) + y;

    if (~isempty(varargin))
        h = plot(xx, yy, varargin{1});
    else    
        h = plot(xx, yy);
    end