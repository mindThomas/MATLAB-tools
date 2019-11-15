function immarker(x, y, size, args, width)
    % Insert a marker on am image
    type = '+';
    color = 'r';    
            
    if (isempty(size))
        size = 10;
    end
    if (isempty(width))
        width = 2;
    end
    if (ischar(args))
        for (i = 1:length(args))
            if (args(i) == '+' || ...
                args(i) == 'x' || ...
                args(i) == 'o' )
                type = args(i);
            end
            if (args(i) == 'r' || ...
                args(i) == 'g' || ...
                args(i) == 'b' || ...
                args(i) == 'c' || ...
                args(i) == 'k' || ...
                args(i) == 'm' || ...
                args(i) == 'y' )
                color = args(i);
            end
        end
    end

    x0 = x - size;
    x1 = x + size;
    y0 = y - size;
    y1 = y + size;    
    if (type == '+')
        line([x0,x1], [y,y], 'Color', color, 'LineWidth', width)
        line([x,x], [y0,y1], 'Color', color, 'LineWidth', width)
    elseif (type == 'x')
        line([x0,x1], [y0,y1], 'Color', color, 'LineWidth', width)
        line([x0,x1], [y1,y0], 'Color', color, 'LineWidth', width)
    elseif (type == 'o')        
        rectangle('Position',[x0,y0, 2*size,2*size],'Curvature',[1 1],'EdgeColor',color,'LineWidth',width)
    end
end