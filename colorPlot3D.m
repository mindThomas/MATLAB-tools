function colorPlot2D(x, y, z, color)
    if (length(x) ~= length(y) || length(x) ~= length(z) || length(x) ~= length(color))
        error('Vector length mismatch');;
    end

    if (size(x, 1) > 1)
        x = x';
    end
    if (size(y, 1) > 1)
        y = y';
    end
    if (size(z, 1) > 1)
        z = z';
    end    
    if (size(color, 1) > 1)
        color = color';
    end    

    if (size(x, 1) ~= 1 || size(y, 1) ~= 1 || size(z, 1) ~= 1 || size(color, 1) ~= 1)
        error('Incorrect input size');
    end

    surface([x;x],[y;y],[z;z],[color;color],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',2);    
end