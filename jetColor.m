function [R,G,B] = jetColor(v, vmin, vmax)
    % White initially
    R = ones(size(v));
    G = ones(size(v));
    B = ones(size(v));
    
    v(find(v < vmin)) = vmin;
    v(find(v > vmax)) = vmax;   
    dv = vmax - vmin;   
    
    for (i = 1:length(v))
    if (v(i) < (vmin + 0.25 * dv))
        R(i) = 0;
        G(i) = 4 * (v(i) - vmin) / dv;
    elseif (v(i) < (vmin + 0.5 * dv))
        R(i) = 0;
        B(i) = 1 + 4 * (vmin + 0.25 * dv - v(i)) / dv;
    elseif (v(i) < (vmin + 0.75 * dv))
        R(i) = 4 * (v(i) - vmin - 0.5 * dv) / dv;
        B(i) = 0;
    else
        G(i) = 1 + 4 * (vmin + 0.75 * dv - v(i)) / dv;
        B(i) = 0;
    end
        
end