function intersection = LineIntersection(line1_start, line1_end, line2_start, line2_end)
    x1 = line1_start(1);
    y1 = line1_start(2);
    x2 = line1_end(1);
    y2 = line1_end(2);
    x3 = line2_start(1);
    y3 = line2_start(2);
    x4 = line2_end(1);
    y4 = line2_end(2);
    
    D = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4);
    intersection = [nan;nan];
    
    if D ~= 0
        intersection = [
            (x1*y2-y1*x2)*(x3-x4) - (x1-x2)*(x3*y4-y3*x4);
            (x1*y2-y1*x2)*(y3-y4) - (y1-y2)*(x3*y4-y3*x4);
        ] / D;
    end
end